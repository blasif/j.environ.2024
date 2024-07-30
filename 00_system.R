
tmp_dir <- getwd()

suppressMessages(require(cocons))
suppressMessages(require(maps))
suppressMessages(require(fields))
suppressMessages(require(ncdf4))
suppressMessages(require(elevatr))
suppressMessages(require(sf))
suppressMessages(require(sp))
suppressMessages(require(geodata))

RScript <- 'Rscript'

if (version$major != '4'){
  stop('use R-4.x.y')
}

ncores <- floor((parallel::detectCores() - 1) / 2) + 1

master_mars <- c(4.5, 4.5, 0.5, 0.5)
master_mars_second <- master_mars
master_mars_second[c(2,4)] <- c(0,2)
master_oma <- c(0, 0, 0, 0)
master_width <- 3 * 1.5
master_height <- 3
master_res <- 150