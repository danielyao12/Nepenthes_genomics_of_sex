# bootstrap_test_equalise_sample_sizes.R
# R code to test statistical significance of differences in mean for two samples with vastly different sizes by means of bootstrap and permutation
# Copyright (C) 2017, ETH Zurich, Mathias Scharmann
# 
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
# 		
# If you use this code please cite:
#
# "Scharmann M, Grafe TU, Metali F, Widmer A. (2017) Sex-determination 
# and sex chromosomes are shared across the radiation of dioecious 
# Nepenthes pitcher plants. XXX"
# 	
# contact: mathias.scharmann[-at-]env.ethz.ch or msph52[-at-]gmail.com


bootstrap_test_equalise_sample_sizes <- function(x,y,n_boostraps) {

# tests for difference in means, bootstrap resampled data have IDENTICAL length in both of the vectors to be compared!
# Null hypothesis: mean(x) - mean(y) = 0

resampling_size = min( length(x), length(y) )

collector <- vector(mode="numeric", length=n_boostraps)
for (i in seq(n_boostraps)){
  
  # bootstrap observations
  idxes = sample(seq( length(x) ), resampling_size, replace = TRUE)
  a = x[idxes]
  idxes = sample(seq( length(y) ), resampling_size, replace = TRUE)
  b = y[idxes]
  
  # now permute
  c = sample( c(a,b) ) # random order of bootstrapped observations mixed at equal weight from both groups
  p_x = c[ 1:resampling_size]
  p_y = c[ (resampling_size+1):length(c) ]
   
  # get difference in means:
  diff_in_means = mean(p_x) - mean(p_y)
  
  collector[i] <- abs( diff_in_means )
  
}

obs_diff = mean(x) - mean(y)
two_sided_excess = collector[ collector >= abs(obs_diff)]

p = max( (length(two_sided_excess)/n_boostraps ) , (1/n_boostraps) )

results = c( mean(collector), sd(collector), obs_diff, p )
names(results) = c("mean_of_mean_difference_of_resampled_data", "SD", "observed_difference_in_means", "p_value")
return( results)

}

# example:
bootstrap_test_equalise_sample_sizes(numeric_vector_1,numeric_vector_2,100000)
