ls.obj <- {as.data.table(sapply(ls(),
     function(x){format(object.size(get(x)),
        nsmall=3,digits=3,unit="Mb")}),keep.rownames=TRUE)[,
        c("mem","unit") := tstrsplit(V2, " ", fixed=TRUE)][,
        setnames(.SD,"V1","obj")][,.(obj,mem=as.numeric(mem),unit)][order(-mem)]}

ls.obj


