library(juniper0)
library('getopt')


spec = matrix(c(
  'input'     , 'i', 1, "character",
  'outprefix' , 'o', 1, "character",
  'mu'        , 'm', 1, "double",
  'serial'    , 's', 1, "double",
  'popsize'   , 'n', 1, "double",
  'fixmu'     , 'M', 0, "logical"
), byrow=TRUE, ncol=4)

opt = getopt(spec)

if (is.null(opt$input)) opt$input <- "./input_data/"
if (is.null(opt$input)) opt$input <- "./juniper"
#if (is.null(opt$mu)) opt$mu <- NA
if (is.null(opt$fixmu)) opt$fixmu <- FALSE
if (is.null(opt$serial)) opt$serial <- 5

init <- initialize(indir  =opt$input,
                  a_g     =opt$serial,
                  init_mu =opt$mu,
                  fixed_mu=opt$fixmu,
                  init_pi =0.05,
                  filters=list(af=0.05, dp=100, sb=10))

outbreak <- run_mcmc(init)

summary <- juniper0::summarize(outbreak)

save(summary,file=paste(outprefix,"summary.Rdata",sep="_"))

