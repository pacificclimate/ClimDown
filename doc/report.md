# Rewriting QDM

As part of this project, PCIC put a large effort in to refactoring and rewriting the research software which executes QDM. There were a number of goals for this exercise, including but not limited to:

* to make the software more widely deployable in production
* to increase robustness
* to decrease or eliminate operator intervention
* to generally decrease the computational requirements (both memory and processor time), but also
* to make the computational requirements more configurable (i.e. to increase computational efficiency when greater resources are available)

## Creation of the `ClimDown` R package

Up until this point, the climate downscaling pipeline has been
written as a loosely affiliated collection of R scripts that had to be
manually run in a particular sequence. The data processing pipeline
required multiple operator intervention steps between each script and
a significant amount of intermediate data storage.

We have reorganized and refactored the code base to be installed and
run as a standard R package. This allows anyone on any machine to
install and load the climate downscaling code with two commands:

~~~R
    > install.packages('ClimDown_0.0.1.tar.gz')
    > library(ClimDown)
~~~

The entire package has been released [on
GitHub](https://github.com/pacificclimate/ClimDown) under and [open
source license](http://www.gnu.org/licenses/gpl-3.0.en.html), so that
it's available for other stakeholders to run and to make
contributions.

We have undertaken the process of combining redundant scripts for
different variables into library code which parameterizes variables
and can be run for any arbitrary variable. Parameterizing variables
has eliminated a lot of redundancy and follows the "Don't Repeat
Yourself (DRY)" Principle.

We have eliminated numerous multi-step downscaling pipelines. For
example the Bias Correction Climate Imprint (BCCI) was formerly five
separate scripts/steps which we have combined into one single wrapper
function that can also be executed as a script. Bias Correction
Constructed Analogues (BCCA) was four scripts which have been combined
into one. The Quantile Perturbation Quantile Mapping (QPQM) was two
different scripts which have been combined into one.

Additionally, we have completely refactored the user configuration for
the downscaling process. At the beginning of this project we started
with 500 LOC of hard-coded configuration and queue submission ("PBS")
scripts which would have to be edited/re-executed for every single
invocation. We have replaced this with 50 LOCs of config code that
utilizes R's built-in
[`options()`](https://stat.ethz.ch/R-manual/R-devel/library/base/html/options.html)
method. This provides sensible defaults, but gives users the option of
configuring the package on a per-system or per-run basis with one
optional line of config in the user's home directory (or sourced from
environment variables).

For example, a user on a system for 4GB of available memory and 8
available cores can simply have a one-line `$HOME/.Rprofile` file which
says:

~~~R
options(max.GB=4, mc.cores=8, par.backend='multicore')
~~~

This configuration method provides both power and flexibility to
accommodate a variety of systems, resource levels, and operators.

Overall, at the start of the project we began with twenty-five
different R files and we have refactored this set down to a single
config file, five library files and four executable scripts. We
started with approximately 3000 lines of code (LOCs) and were able to
refactor down to about 1200 LOCs. This simplification makes the
package significantly more accessible to new scientific users.

### `udunits2` and `PCICt`

We eliminated a substantial amount of code from the package by
delegating common abstractions to other libraries specially designed
for these tasks. We delegated all units conversion to the `udunits2`
library. This helped to eliminate code, but also provides some speedup
benefits by performing simple, but expensive vector arithmetic at a
lower level in the software stack.

We rewrote all home-grown time handling code to instead utilize the
[`PCICt` R package](https://cran.r-project.org/package=PCICt). This
package stores dates/times in a proper time structure which:

1. Cuts the memory requirements for all time values by 2-3 times
1. Properly handles all calendar variations in one place
1. Provides compatibility with all other time handling code
1. Delegates specialization to a library which is designed to handle it

## Rewriting numerous algorithms

From the previous section, it may be easy to think that the
investments in the ClimDown code base were solely
*aesthetic*. However, we put a significant effort into rewriting the
core algorithms of BCCA, BCCI and QPQM to achieve greater
computational, memory and storage efficiency. In the following
subsections we describe a few of the more significant successes of
that effort.

### BCCA 47% speedup on BCCA

Bias Correction Constructed Analogues (BCCA) is a particularly
expensive part of the processing chain, requiring upwards of 70 hours
of processing time in some cases. Overall, we were able to realize
approximately a 2x speedup using a number of optimizations. Some of
the more particularly successful methods are described in the
subsections below.

### 14,000x speedup of BCCA grid mapping

The initial step in BCCA involves fine scale grid cells to coarse
scale grid cells and aggregating the values of the fine to the
coarse. For example, say the global grid is 50x25 and the local grid
is 1000x500. For each grid cell in the local grid, BCCA determines
which grid cell in the global grid it corresponds.

The naïve way to think about this problem would be to try and minimize
the distance between L[n] and G[n]. So a simple way to do the search
would be:

~~~
for each Local cell L[i]:
  for each Global cell G[j]:
     compute distance between L[i] and G[j]
  find the minimum distance in the set L[i] * G
  return the index of the minimum
~~~

This seems simple enough, However, the careful observer will note that
the processor has to do a *lot* of extra work. Examining the algorithm
in terms of the size of the input we see that:

~~~
for each Local cell L[i]:                        # Do this L times
  for each Global cell G[j]:                     # Do this L x G times
     compute distance (d) between L[i] and G[j]  # Do this L x G times
  find the minimum distance in the set d[i*j]    # Read G cells L times (cost L x G)
  find the index whose cell matches the minimum  # Read G cells L times (cost L x G)
~~~

Formerly, the [code for
this](https://github.com/pacificclimate/ClimDown/blob/45728050cd5d12ac6703286df4d67f040f7a1cb3/BCCA/BCCA.1.setup_file.R#L77)
looked something like this:

~~~R
obs.lon <- ncvar_get(nc.obs, 'lon')
obs.lat <- ncvar_get(nc.obs, 'lat')
n.lon <- length(obs.lon)
n.lat <- length(obs.lat)

obs.lats <- matrix(obs.lat, nrow=n.lon, ncol=n.lat, byrow=TRUE)
obs.lons <- matrix(obs.lon, nrow=n.lon, ncol=n.lat)
obs.time <- netcdf.calendar(nc.obs)

gcm.lon <- ncvar_get(nc.gcm, 'lon')-360
gcm.lat <- ncvar_get(nc.gcm, 'lat')
gcm.lats <- matrix(gcm.lat, ncol=length(gcm.lat), nrow=length(gcm.lon),
                   byrow=TRUE)
gcm.lons <- matrix(gcm.lon, ncol=length(gcm.lat), nrow=length(gcm.lon))
gcm.lons.lats <- cbind(c(gcm.lons), c(gcm.lats))

# Figure out which GCM grid boxes are associated with each fine-scale grid point
# Confine search to 10 deg. x 10 deg. neighbourhood

dxy <- 10
mdist <- function(x, y)
    apply(abs(sweep(data.matrix(y), 2, data.matrix(x), '-')), 1, sum)
nn <- list()
for (i in seq_along(obs.lons)) {
    if((i %% 500)==0) cat(i, '')
    gcm.lims <- ((gcm.lons.lats[,1] >= (obs.lons[i]-dxy)) &
                 (gcm.lons.lats[,1] <= (obs.lons[i]+dxy))) &
                ((gcm.lons.lats[,2] >= (obs.lats[i]-dxy)) &
                 (gcm.lons.lats[,2] <= (obs.lats[i]+dxy)))
    gcm.lims <- which(gcm.lims)
    nn.min <- which.min(mdist(c(obs.lons[i], obs.lats[i]),
                        gcm.lons.lats[gcm.lims,]))
    nn[[i]] <- gcm.lims[nn.min]
}
nn <- unlist(nn)
~~~

This seemed like a simple algorithm: just compute the cell distances
and then find the minimum. However, it was also problematic because as
the size of the number of local cells grows, our cost of computation
grows by its product with the number of global grid cells. For
Canadian ANUSPLIN data, there are 1068 x 510 cells (for a total of
544,680) and let's say that our GCM has 50 x 25 cells (for a total of
1,250 cells). So the cost of the inner loop in "some computational
unit" is:

$$ (c_0 \cdot L \times G) + (c_1 \cdot  L \times G) + (c_2 \cdot L \times G) $$

where the $c$ terms are constants that correspond to the cost of
computing a distance between two points, finding the minimum point,
and finding an array index. The constant terms are not important,
since they are not affected by the size of the input. So we can just
group them together and call the cost;

$$ (c \cdot L \times G) $$

So for this set of input, our cost is $$ 544,680 \times 1,250 =
680,850,000 $$ or approximately 680 million "computational
units". Running the naïve implementation took about 1668 seconds which
is a little less than half an hour.

~~~R
> source('BCCA/naive.implementation.R')
500 1000 1500 2000 2500 3000 ... 543000 543500 544000 544500 [1] "Elapsed Time"
    user   system  elapsed 
1668.868    8.926 1681.728 
~~~

In order to improve upon this algorithm, one would want to take
advantage of the shared structure between the two grids that are being
compared. For example the latitudes and longitudes in both the coarse
and the fine grid are in sorted order. In order to search for a
number, one doesn't have to look at every single number. We can use a
bisect algorithm that examines the point in the middle and then
decides which half of the array to search. Then searching the full
space only costs you the log (base 2) of the search space.

The other major structure that we hadn't been taken advantage of is
the fact that the latitudes repeat themselves in the $x$ dimension and
the longitudes repeat themselves in the $y$ dimension. So instead of
doing an operation $x \times y$ times, we can do it $x + y$
times. That's a *huge* optimization.

Pseudo-code for such a search would look like this:

~~~
For each local[x]:
    bisect_search(local[x], Global[x])

For each local[y]:
    bisect_search(local[y], Global[y])

return a 2d grid of the search results for each dimension
~~~

In code (now [found in the package](https://github.com/pacificclimate/ClimDown/blob/master/R/bisect.R)):

~~~R
## Perform a binary search on the *sorted* vector v
## Return the array index of the element closest to x
find.nearest <- function(x, v) {
    if (length(v) == 1) {
        return(1)
    }
    if (length(v) == 2) {
        return(which.min(abs(v - x)))
    }
    mid <- ceiling(length(v) / 2)
    if (x == v[mid]) {
        return(mid)
    } else if (x < v[mid]) {
        return(find.nearest(x, v[1:mid]))
    }
    else {
        return((mid - 1) + find.nearest(x, v[mid:length(v)]))
    }
}

regrid.one.dim <- function(coarse.points, fine.points) {
    return(sapply(fine.points, find.nearest, coarse.points))
}

## Take a fine scale (e.g. ANUSPLIN) grid of latitudes and longitudes
## and find the indicies that correspond to a coarse scale (e.g. a GCM) grid
## Since the search is essentially a minimizing distance in 2 dimensions
## We can actually search independently in each dimensions separately (which
## is a huge optimization, making the run time x + y instead of x * y) and
## then reconstruct the indices to create a full grid
regrid.coarse.to.fine <- function(coarse.lats, coarse.lons, fine.lats, fine.lons) {
    xi <- regrid.one.dim(gcm.lon, obs.lon)
    yi <- regrid.one.dim(gcm.lat, obs.lat)
    ## Two dimensional grid of indices
    xi <- matrix(xi, ncol=length(fine.lats), nrow=length(fine.lons), byrow=F)
    yi <- matrix(yi, ncol=length(fine.lats), nrow=length(fine.lons), byrow=T)
    return(list(xi=xi, yi=yi))
}

~~~

The cost for every bisection search is the log of the input size. Our
input size is divided into X and Y space this time, so this can be
represented by $G_x, G_y, L_x$, and $L_y$ for Global, Local, X and Y.

$$
cost =  L_x \times log_2 G_x + L_y \times log_2 G_y + L_x \times L_y
$$

Substituting the input numbers into this equation gives us a cost
estimate of 553,076 (compared to 680 million). We were able to observe
this speedup in the run time.

~~~R
> ptm <- proc.time()
> rv <- regrid.coarse.to.fine(gcm.lat, gcm.lon, obs.lat, obs.lon)
> print('Elapsed Time')
> print(proc.time() - ptm)
[1] "Elapsed Time"
   user  system elapsed 
  0.117   0.000   0.117 
> str(rv)
List of 2
 $ xi: num [1:1068, 1:510] 15 15 15 15 15 15 15 15 15 15 ...
 $ yi: num [1:1068, 1:510] 13 13 13 13 13 13 13 13 13 13 ...
~~~

What had previous taken us nearly half an hour to compute, we were
able to reduce to slightly over $\frac{1}{10}$ of a second, achieving
a 14,000x speedup.

~~~R
> 1668.868 / .117
[1] 14263.83
~~~

It may be noted that this step didn't necessarily consist a large
percentage of QDM's overall run time. However, 30 minutes was
significant enough that this step's output had to be saved to disk and
manually checked by an operator before proceeding. By making this
step essentially free to compute, we were able to integrate it
directly into the BCCA pipeline, eliminate operator intervention and
make the whole process much more efficient.

### Selecting Analogues

Selecting analogues involves choosing the nearest $k$ timesteps
(generally $k=30$, but this is globally configurable by saying
`options(n.analogues=k)`), or the $k$ timesteps with the smallest
squared difference to the current timestep.

A naïve implementation of selecting the smallest differences is to
take the compute the difference of every timestep, sort the
differences and take the first $k$ values in the sorted list.

A potential downside of this method is that sorting an entire list
requires a computational complexity of $O(n \cdot log_2 n)$ and one wastes time
on ordering the elements of the list in which we have no interest
(i.e. all those that are greater than the $kth$ element).

Instead, one can utilize the
[*quickselect*](https://en.wikipedia.org/wiki/Quickselect) algorithm
which has an average performance of $O(n)$. `quickselect` selects the
$k^{th}$ largest element in a series, then one can run a final pass
over the series to select all elements that are smaller than the
$k^{th}$ largest.

For large input sizes (i.e. for long timeseries) `quickselect` performs
much, much better than the naïve implementation of simply
sorting. However, sorting actually has a relatively low constant
cost. So for smaller input sizes, sorting actually wins.

We have experimented with input sizes and determined 1 million input
points to be the approximate point at which the curves cross. We have
codified this in
[`quickselect.R`](https://github.com/pacificclimate/ClimDown/blob/master/R/quickselect.R). Since
our current input sets run at about 55,000 timesteps, the sorting
method will be used for the foreseeable future. However, if the code
ends up being run on sub-daily output, or multi-century runs, the
package will automatically do the right thing and select the faster
algorithm.

While this effort has not necessarily increased the *speed* of this
section of the code, it has significantly increased the *scalability* of
the code for longer or higher-resolution timeseries.

### 12,000x reduction of BCCA output space

As mentioned in the previous section, the BCCA process involves
selecting the top 30 timesteps with the smallest squared differences
to a given timestep. Each of these 30 steps are then assigned a
fractional weights, and the new timestep calculated as a weighted sum.

In the original implementation of BCCA, the final output consisted of
the full high spatial resolution output with the constructed analogues
applied. This required $1068 \times 510 \text{ cells} \times 54750
\text{ timesteps} \times 4 \text{ bytes}$ (approximately 56 GB) of
intermediate storage for each model.

An optimization which can be performed is to defer the application of
the analogues until later in the processing chain. Doing so means that
BCCA only needs to output $k$ timestep indices (i.e. the analogues)
plus $k$ weights for each timestep. The indices can be stored as 2
byte short integers while the weights can be stored as a 4 byte float.

Using this method, rather than the storage costs being:

$$ 4 \text{ bytes} \times c \text{ cells} \times t \text{ timesteps} $$

the storage costs end up being:

$$ 6 \text{ bytes} \times t \text{ timesteps} \times k \text{ analogues} $$

In our case, $c$ is extremely large (544,680), and $k$ is small (30)
so we end up reducing the storage requirements by $544,680 / 30 /
\frac{3}{2} = 12,104$ times. In the end, this scheme only uses about 5
MB, which is small enough to easily store in RAM and simply pass along
to the next step in the processing chain. This saves us a significant
amount of time I/O time reading and writing to disk (on the order of
20 minutes for every variable on every model) and a significant amount
of intermediary storage (56 GB on every variable on every model).

### 5x speedup of QPQM

The Quantile perturbation quantile mapping (QPQM) algorithm is one
stage of the algorithm that takes a significant amount of time. We
performed aggressive profiling and rewrote everything in the code, but
the core computation. We were able to achieve a 4-5 times speedup from
the code that existed at the beginning of this project.

The improvements that we made include the following:

* We rewrote the chunking mechanism to automatically fill the amount
  of available RAM, increasing the I/O efficiency

* We rewrote the time handling code to use proper date times, reducing
  the information in memory

* We converted the date handling code to use R `factors`, rather than
  using repeated expensive comparisons

* We moved all temporal factor calculation *outside* of the main
  spatial loop to avoid recomputing the same thing (this alone cut 1/3
  of the run time)

* We [eliminated certain
  scenarios](https://github.com/pacificclimate/ClimDown/compare/22c369d697bba9f240a4...6ef0e9c4fa2b3f87fa#diff-722dd0bb08ce9a712b4d1386a001d5d4L190)
  where QPQM was calculating three times the volume of results
  necessary and then throwing two thirds of them away.

Finally, after refactoring QPQM to be five times faster, we were able
to easily parallelize the remaining CPU bound loop to achieve greater
speedups when a parallel compute architecture (either MPI or
multicore) is available.