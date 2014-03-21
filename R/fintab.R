fintab <-
function(gen1) {
        k = nrow(gen1)
        k3 = k/2
        k1 = 2 * ncol(gen1)
        k2 = k1/2
        t5 = vector("list", k2)
        t1 = gen1[1:(k/2), ]
        t1 = data.frame(t1)
        t2 = gen1[((k/2) + 1):k, ]
        t2 = data.frame(t2)
        t3 = cbind(t1, t2)
        t3 = data.frame(t3)
        for (i in 1:k2) t5[[i]] = cbind(t3[, i], t3[, i + k2])
        t6 = vector("list", k2)
        for (i in 1:k2) t6[[i]] = matrix(0, nrow = k3, ncol = 2)
        for (j in 1:k3) {
            for (i in 1:k2) t6[[i]][j, ] = t(as.data.frame(range((t5[[i]])[j, 
                ])))
        }
        t7 = t6[[1]]
        {
            for (s in 2:k2) t7 = cbind(t7, t6[[s]])
        }
        t7 = data.frame(t7)
        return(t7)
    }
