vonMises <- R6Class("vonMises",
                    public = list( 
                      initialize=function(mean,kappa){
                        if(!missing(mean)) private$mean <- mean
                        if(!missing(kappa)) private$kappa <- kappa
                      },
                      set_mean = function(val){
                        private$mean <- val
                      },
                      set_kappa=function(val){
                        private$kappa<-val
                      },
                      get_mean = function()return (private$mean), 
                      get_kappa=function() return (private$kappa),
                      density = function(values,loglik) return(dvonmises(values, mu=private$mean, kappa=private$kappa, log=loglik))
                    ),
                    private = list(
                      mean=0,
                      kappa=1
                    )
)
