library(akpvdata)
data(akpvpolys)

library(JHop)
data(dHOfal)
polyids = unique(dHOfal$polyid)
plot(akpvpolys, xlim = c(950125.1, 1289952.4),
	ylim = c(1000471.7, 1178765.9))
plot(akpvpolys[akpvpolys$polyid %in% polyids[1],], 
	col = 'red', 	add = TRUE)
plot(akpvpolys[akpvpolys$polyid %in% polyids[2],], 
	col = 'blue', 	add = TRUE)
plot(akpvpolys[akpvpolys$polyid %in% polyids[3],], 
	col = 'green', 	add = TRUE)
plot(akpvpolys[akpvpolys$polyid %in% polyids[4],], 
	col = 'orange', 	add = TRUE)

# locator()

polyids[1]
sum(dHOfal$polyid == polyids[1])
sum(dHOfal$polyid == polyids[2])
sum(dHOfal$polyid == polyids[3])
sum(dHOfal$polyid == polyids[4])
