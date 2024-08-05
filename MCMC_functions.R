#### Please do not run this file directly as it will not work.
#### Please refer to file 00_run.R

sampleDistGaussian<- function(mean, sd, N=1, min=0, max=1){
	repeat{
		point= rnorm(N,mean=mean,sd=sd)
		if(point>=min & point<=max) break;
	}
	return(point)
}