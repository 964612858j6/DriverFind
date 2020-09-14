#Name: DriverFind
#Author: Mingyuan Luan, Rongfeng Huang
#Maintainer: Mingyuan Luan <964612858@qq.com>

if(!require("data.table")) install.packages("data.table",update = F,ask = F)
if(!require("parallel")) install.packages("parallel",update = F,ask = F)
if(!require("splitstackshape")) install.packages("splitstackshape",update = F,ask = F)
if(!require("GSVA")) install.packages("GSVA",update = F,ask = F)
if(!require("PMCMRplus")) install.packages("PMCMRplus",update = F,ask = F)
if(!require("ggplot2")) install.packages("ggplot2",update = F,ask = F)
	
DriverFind=function(omics=NULL, category=NULL, cores=NULL, gsva_method="gsva", pathway_list=NULL){

	if(!require("data.table")){
		stop("R package data.table should be installed.")
	}
	if(!require("parallel")){
		stop("R package parallel should be installed.")
	}
	if(!require("splitstackshape")){
		stop("R package splitstackshape should be installed.")
	}
	if(!require("GSVA")){
		stop("R package GSVA should be installed.")
	}
	if(!require("PMCMRplus")){
		stop("R package PMCMRplus should be installed.")
	}
	if(!require("ggplot2")){
		stop("R package ggplot2 should be installed.")
	}

	if (is.null(omics)){
		stop("Omics data should be provided.")
	}
	if (is.null(category)){
		stop("The category of omics data should be provided.")
	}
	if (length(category)!=length(omics)){
		stop("The category of all omics data should be provided.")
	}
	if (any(is.na(match(category, c("mrna", "cnv", "cna", "cnd", "met", "mut"))))){
		err=category[which(is.na(match(category, c("mrna", "cna", "cnd", "met", "mut"))))]
		ret="Please check the spelling of category: "
		ret=paste0(ret, glue(err, sep=", "))
		stop(ret)
	}
	if (length(which(category=="mrna"))==0){
		stop("mRNA data should be provided.")
	}
	if (length(category)<2){
		stop("Less than two omics data can not run this algorithm.")
	}
	if (is.null(cores)){
		cores=detectCores()
	}
	if (is.null(names(pathway_list))){
		stop("pathway names in pathway_list should be provided.")
	}


	####Support function 01  path_alter_find
	path_alter_find<<-function(path_name){
		path_name=as.character(as.matrix(path_name))
		###通路对应的基因名单
		genelist=out[[which(name==path_name)]]
		
		alter=alist()
		for (i in 1:length(omics_new)){
			tag=na.omit(match(genelist, rownames(omics_new[[i]])))
			if (length(tag)>0){
				alter[[i]]=t(omics_new[[i]][tag, ])
				if (nrow(alter[[i]])==1){
					alter[[i]]=t(alter[[i]])
					colnames(alter[[i]])=rownames(omics_new[[i]])[tag]
				}
			}
			else{
				alter[[i]]=NULL
			}
		}
		return(alter)
	}
	####Support function 02  score_compar
	score_compar<<-function(alter, score, type, split_by="!!@@"){
		
		######Support function
		###NO_01 glue
		glue<-function(x, sep=NULL){
			if (is.null(sep)){
				for (i in 1:length(x)){
					if (i==1){
						y=x[1]
					}
					else{
						y=paste0(y, x[i])
					}
				}
			}
			
			if (!is.null(sep)){
				for (i in 1:length(x)){
					if (i==1){
						y=x[1]
					}
					else{
						y=paste(y, x[i], sep=sep)
					}
				}
			}
			
			return(y)
		}
		
		###NO_02 fdr_filter
		fdr_filter<-function(dat, fdr_col=NULL, criterion=NULL){
			as_character=function(x){return(as.character(as.matrix(x)))}
			as_numeric=function(x){return(as.numeric(as.matrix(x)))}

			dat=as.data.frame(dat)
			
			if (!is.null(fdr_col)){
				if (class(try(fdr_col%%1==0, silent = TRUE))=="try-error"){
					fdr_col=which(colnames(dat)==fdr_col)
				}
				
				p=as_numeric(dat[, fdr_col])
				p=cbind(p, 1:length(p))
				colnames(p)=c("V1", "V2")
				p=p[order(p[,1]), ]
				q=p.adjust(as_numeric(p[,1]), method="fdr")
				p[,1]=q
				p=as_numeric(p[order(as_numeric(p[, 2])), 1])
				dat[, fdr_col]=p
				if (!is.null(criterion)){
					dat=dat[which(p<criterion), ]
				}
				return(dat)
			}
			else{
				p=as_numeric(dat)
				p=cbind(p, 1:length(p))
				colnames(p)=c("V1", "V2")
				p=p[order(p[,1]), ]
				q=p.adjust(as_numeric(p[,1]), method="fdr")
				p[,1]=q
				
				if (!is.null(criterion)){
					p=p[order(as_numeric(p[, 2])), ]
					p=p[which(p[,1]<criterion), ]
					return(p)
				}
				else{
					p=as_numeric(p[order(as_numeric(p[, 2])), 1])
					return(p)
				}
			}
		}
		
		split_by <- split_by
		type <- type
		
		score=as.data.frame(score)
		rownames(score)=as.character(as.matrix(score[,1]))
		alter[[length(alter)+1]]=score
		
		###find which is null
		for (i in 1:length(alter)){
			if (i==1){
				nu=is.null(alter[[i]])
			}
			else{
				nu=c(nu, is.null(alter[[i]]))
			}
		}
		
		for (i in 1:length(which(!nu))){
			if (i==1){
				tag=rownames(alter[[which(!nu)[i]]])
			}
			else{
				tag=intersect(tag, rownames(alter[[which(!nu)[i]]]))
			}
		}
		
		
		coln=alist()
		for (i in which(!nu)){
			coln[[i]]=colnames(alter[[i]])
		}
		for (i in which(!nu)){
			
			alter[[i]]=alter[[i]][match(tag, rownames(alter[[i]])), ]
			if (all(is.null(dim(alter[[i]])), length(alter[[i]])>0)){
				alter[[i]]=as.data.frame(alter[[i]])
			}
			colnames(alter[[i]])=coln[[i]]
		}
		
		if (any(nu)){
			alter[[which(nu)]]="NULL"
		}
		
		sc<-alter[[length(alter)]]
		sc=as.data.frame(sc); sc[,2]=as.numeric(as.matrix(sc[,2]))
		sc<<-sc
		
		alt<<-alter[1:(length(alter)-1)]
		
		##summary the condition for calculation of mean
		tong=function(dat){
			return(names(table(as.character(as.matrix(apply(dat, 2, as.character))))))
		}
		
		all=alist()
		for (i in 1:length(alt)){
			if (class(try(tong(alt[[i]]), silent = TRUE))=="try-error"){
				all[[i]]=NULL
			}
			else{
				all[[i]]=tong(alt[[i]])
			}
		}
		for (i in 1:length(all)){
			if (is.null(all[[i]])){
				if (i==1){
					alls=NA
				}
				else{
					alls=c(alls, NA)
				}
				next
			}
			if (i==1){
				alls=as.character(as.matrix(all[[i]]))
			}
			else{
				alls=c(alls, as.character(as.matrix(all[[i]])))
			}
		}
		rm(all); gc()
		ss=as.character(as.matrix(na.omit(alls)))
		rm(alls); gc()
		alls<<-names(table(ss))
		
		###filter genes with less than 5% altertions.
		fun=function(dat){
			sss<<-dat
			core=function(i){
				x=as.character(as.matrix(sss[,i]))
				return(length(which(x=="1")))
			}
			return(sapply(1:ncol(sss), core))
		}
		filt=alist()
		nu=nu[1:length(alt)]
		cutof=0.05*nrow(alt[[which(!nu)[1]]])
		
		for (i in which(!nu)){
			filt[[i]]=fun(alt[[i]])
			len=which(filt[[i]]>cutof)
			if (length(len)>0){
				if (length(len)==1){
					nam=colnames(alt[[i]])[len]
					alt[[i]]=as.data.frame(alt[[i]][, len])
					colnames(alt[[i]])=nam
				}
				else{
					alt[[i]]=alt[[i]][, len]
				}
			}
			else{
				alt[[i]]="NULL"
			}
		}
		
		
		################Only three groups calculated p_s, both 2 groups and three groups calculated the mean values in each status.
		out=alist()
		genename=alist()
		sel=vector(length=length(alt))
		for (j in 1:length(alt)){
				
			mat<-alt[[j]]
			tt<-type[j]

			
			if (any(mat=="NULL")){
				sout=alist()
				sout[[1]]=glue(c(tt, NA, NA), sep=split_by)
				
				###p_s
				p_s=NA
				sout[[2]]=p_s
				
				###me
				me=glue(rep(NA, length(alls)), sep=split_by)
				sout[[3]]=me
				
				out[[j]]=sout
				
				genename[[j]]=NA
				sel[j]="yes"
				next
			}
			else{
				sel[j]="no"
				inter_outcome=alist()
				for (i in 1:ncol(mat)){
					tag=as.character(as.matrix(mat[, i]))
					shan=table(tag)
					if (length(shan)<2){
						com=glue(c(tt, NA, NA), sep=split_by)
						p_s=NA
						me=glue(rep(NA, length(alls)), sep=split_by)
						
						s_outcome=alist()
						s_outcome[[1]]=com; s_outcome[[2]]=p_s; s_outcome[[3]]=me
						inter_outcome[[i]]=s_outcome
					}
					if (length(shan)==2){
						x1=as.numeric(as.matrix(sc[which(tag==names(shan)[1]),2]))
						x2=as.numeric(as.matrix(sc[which(tag==names(shan)[2]),2]))
						com=glue(c(tt, "wilcox", wilcox.test(x2, x1, exact=FALSE)$p.value), sep=split_by)
						
						###p_s
						p_s=NA
						
						###me
						##There is no need to worry about the order, x1 and x2 are ordered by the order of shan
						for (w in 1:length(shan)){
							if (w==1){
								me=mean(x1)
							}
							else{
								me=c(me, mean(x2))
							}
						}
						names(me)=names(shan)
						me=me[match(alls, names(me))]
						me=glue(me, sep=split_by)
						
						s_outcome=alist()
						s_outcome[[1]]=com; s_outcome[[2]]=p_s; s_outcome[[3]]=me
						inter_outcome[[i]]=s_outcome
					}
					if (length(shan)>2){
						sc_com=as.numeric(as.matrix(sc[, 2]))
						group=as.factor(tag)
						com=glue(c(tt, "kruskal", kruskal.test(sc_com~group)$p.value), sep=split_by)
						
						###p_s
						p_s=PMCMRplus::kwAllPairsNemenyiTest(sc_com ~ group)
						
						###me
						for (w in 1:length(shan)){
							if (w==1){
								me=mean(sc_com[which(tag==names(shan)[1])])
							}
							else{
								me=c(me, mean(sc_com[which(tag==names(shan)[w])]))
							}
						}
						names(me)=names(shan)
						me=me[match(alls, names(me))]
						me=glue(me, sep=split_by)
						
						s_outcome=alist()
						s_outcome[[1]]=com; s_outcome[[2]]=p_s; s_outcome[[3]]=me
						inter_outcome[[i]]=s_outcome
					}
				}
				out[[j]]=inter_outcome
				genename[[j]]=colnames(mat)
			}
		}
		
		###merge the outcome
		chai=function(dat){
			dat=data.table::as.data.table(dat); colnames(dat)="V1"
			dat=splitstackshape::cSplit(dat, "V1", split_by)	
			return(dat)
		}
		
		if (length(table(sel))==1){
			sel=NULL
		}
		if (length(table(sel))==2){
			sel=which(sel=="yes")
		}
		
		com=alist()
		p_s=alist()
		me=alist()
		
		if (is.null(sel)){
			for (i in 1:length(out)){
				sout=out[[i]]
				scom=alist()
				sp_s=alist()
				sme=alist()
				for (j in 1:length(sout)){
					scom[[j]]=sout[[j]][[1]]
					sp_s[[j]]=sout[[j]][[2]]
					sme[[j]]=sout[[j]][[3]]
				}
				com[[i]]=chai(as.character(scom))
				me[[i]]=chai(as.character(sme))
				if (i==1){
					p_s=sp_s
				}
				else{
					p_s=c(p_s, sp_s)
				}
			}
		}
		else{
			for (i in 1:length(out)){
				sout=out[[i]]
				scom=alist()
				sp_s=alist()
				sme=alist()
				if (length(intersect(i, sel))==0){
					for (j in 1:length(sout)){
						scom[[j]]=sout[[j]][[1]]
						sp_s[[j]]=sout[[j]][[2]]
						sme[[j]]=sout[[j]][[3]]
					}
				}
				else{
					for (j in 1:(length(sout)/3)){
						scom[[j]]=sout[[1]]
						sp_s[[j]]=sout[[2]]
						sme[[j]]=sout[[3]]
					}
				}
				com[[i]]=chai(as.character(scom))
				me[[i]]=chai(as.character(sme))
				if (i==1){
					p_s=sp_s
				}
				else{
					p_s=c(p_s, sp_s)
				}
			}
		}
		
		if (is.null(sel)){
			for (j in 1:length(com)){
				if (j==1){
					if (nrow(com[[j]])==1){
						scom=com[[j]]
					}
					else{
						scom=fdr_filter(dat=com[[j]], fdr_col=3, criterion=NULL)
					}
					sme=me[[j]]
					g=as.character(genename[[j]])
				}
				else{
					if (nrow(com[[j]])==1){
						scom=rbind(scom, com[[j]])
					}
					else{
						scom=rbind(scom, fdr_filter(dat=com[[j]], fdr_col=3, criterion=NULL))
					}
					
					sme=rbind(sme, me[[j]])
					g=c(g, as.character(genename[[j]]))
				}
			}
		}
		else{
			for (j in 1:length(out)){
				if (length(intersect(j, sel))>0){
					if (j==1){
						scom=com[[j]]
						sme=me[[j]]
						g=as.character(genename[[j]])
					}
					else{
						scom=rbind(scom, com[[j]])
						sme=rbind(sme, me[[j]])
						g=c(g, as.character(genename[[j]]))
					}
					next
				}
				else{
				if (j==1){
					if (nrow(com[[j]])==1){
						scom=com[[j]]
					}
					else{
						scom=fdr_filter(dat=com[[j]], fdr_col=3, criterion=NULL)
						dat=com[[j]]; fdr_col=3; criterion=NULL
					}
					sme=me[[j]]
					g=as.character(genename[[j]])
				}
				else{
					if (nrow(com[[j]])==1){
						scom=rbind(scom, com[[j]])
					}
					else{
						scom=rbind(scom, fdr_filter(dat=com[[j]], fdr_col=3, criterion=NULL))
					}
					
					sme=rbind(sme, me[[j]])
					g=c(g, as.character(genename[[j]]))
				}
				}
			}
		}
		
		colnames(sme)=alls
		scom=as.data.frame(cbind(g, scom))
		colnames(scom)=c("gene", "category", "method", "Q_value")
		
		out=alist()
		out[[1]]=scom; out[[2]]=p_s; out[[3]]=sme

		return(out)
		rm(sc); rm(alt); rm(mat); rm(tt); rm(cores); gc()
	}
	####Support function 03  path_driver_find
	path_driver_find<<-function(dat, q_cutoff=0.01, sel_pair, ps_cutoff=0.01, split_by="--", return_mean=FALSE){
		com=dat
		com_sel_by=which(as.numeric(as.matrix(com[[1]]$Q_value))<q_cutoff)
		
		###Support function
		glue=function(x, sep=NULL){
			if (is.null(sep)){
				for (i in 1:length(x)){
					if (i==1){
						y=x[1]
					}
					else{
						y=paste0(y, x[i])
					}
				}
			}
			
			if (!is.null(sep)){
				for (i in 1:length(x)){
					if (i==1){
						y=x[1]
					}
					else{
						y=paste(y, x[i], sep=sep)
					}
				}
			}
			
			return(y)
		}

		if (length(com_sel_by)==0){
			return(c(NA, NA, NA))
		}
		else{
			com_sel<<-alist()
			com_sel[[1]]=com[[1]][com_sel_by, ]
			com_sel[[2]]=com[[2]][com_sel_by]
			com_sel[[3]]=com[[3]][com_sel_by, ]

			genename=as.character(as.matrix(com_sel[[1]]$gene))
			category=as.character(as.matrix(com_sel[[1]]$category))
			
			fun=function(i){
				cate=as.character(as.matrix(com_sel[[1]][i, 2]))
				method=as.character(as.matrix(com_sel[[1]][i, 3]))
				#com_sel[[1]][i, ]
				
				
				if (method=="kruskal"){
					p_s<<-com_sel[[2]][[i]]$p.value
					me<<-as.numeric(as.matrix(com_sel[[3]][i, ]))
					names(me)=names(com_sel[[3]])
					
					find_name=function(x1, x2){###比较p_s中的显著性，换做更多组数也适用，无需更改
						x1=as.character(as.matrix(x1))
						x2=as.character(as.matrix(x2))
						loc=alist()
						loc[[1]]=c(which(rownames(p_s)==x2), which(colnames(p_s)==x1))
						loc[[2]]=c(which(rownames(p_s)==x1), which(colnames(p_s)==x2))
						
						tran=function(x){
							if (length(na.omit(x))==2){
								s=p_s[x[1], x[2]]
								if (is.na(s)){
									return(NA)
								}
								else{
									return(s)
								}
							}
							else{
								return(NA)
							}
						}
						loc=as.character(lapply(loc, tran))
						if (length(loc)!=length(which(loc=="NA"))){
							return(loc[which(loc!="NA")])
						}
						else{
							return(NA)
						}
					}
					
					sel_pair=data.table::as.data.table(sel_pair); colnames(sel_pair)="V1"
					sel_pair=splitstackshape::cSplit(sel_pair, "V1", split_by)
					sel_pair=as.data.frame(sel_pair)

					pair=alist()
					for (i in 1:nrow(sel_pair)){
						pair[[i]]=as.numeric(as.matrix(find_name(sel_pair[i,1], sel_pair[i,2])))
					}
					pair=as.numeric(pair)
					name_pair=paste(sel_pair[,1], sel_pair[,2], sep=split_by)
					names(pair)=name_pair
					
					ptag=as.data.frame(sel_pair)
					
					out=alist()
					if (return_mean){
						for (i in 1:nrow(ptag)){
							mename=as.character(as.matrix(ptag[i, ]))
							if(me[which(names(me)==mename[1])]< me[which(names(me)==mename[2])]){
								out[[i]]=glue(c("down", me[which(names(me)==mename[1])], me[which(names(me)==mename[2])]), sep="__")
							}
							if(me[which(names(me)==mename[1])]> me[which(names(me)==mename[2])]){
								out[[i]]=glue(c("up", me[which(names(me)==mename[1])], me[which(names(me)==mename[2])]), sep="__")
							}
							if(me[which(names(me)==mename[1])]== me[which(names(me)==mename[2])]){
								out[[i]]=glue(c("maintain", me[which(names(me)==mename[1])], me[which(names(me)==mename[2])]), sep="__")
							}
						}
					}
					else{
					for (i in 1:nrow(ptag)){
						mename=as.character(as.matrix(ptag[i, ]))
							if(me[which(names(me)==mename[1])]< me[which(names(me)==mename[2])]){
								out[[i]]="down"
							}
							if(me[which(names(me)==mename[1])]> me[which(names(me)==mename[2])]){
								out[[i]]="up"
							}
							if(me[which(names(me)==mename[1])]== me[which(names(me)==mename[2])]){
								out[[i]]="maintain"
							}
						}
					}

					out=as.character(out)
					names(out)=names(pair)
					out[which(pair>=ps_cutoff)]="-"

					return(out)
				}
				if (method=="wilcox"){
					me<<-as.numeric(as.matrix(com_sel[[3]][i, ]))
					names(me)=names(com_sel[[3]])
					
					tt=c(which(sel_pair==paste(1,0,sep=split_by)), which(sel_pair==paste(0,1,sep=split_by)))
					out=rep("-", length(sel_pair))
					
					if (return_mean){
						if(me[which(names(me)==1)]> me[which(names(me)==0)]){
							out[tt]=glue(c("up", me[which(names(me)==1)], me[which(names(me)==0)]), sep="__")
						}
						if(me[which(names(me)==1)]< me[which(names(me)==0)]){
							out[tt]=glue(c("down", me[which(names(me)==1)], me[which(names(me)==0)]), sep="__")
						}
						if(me[which(names(me)==1)]== me[which(names(me)==0)]){
							out[tt]=glue(c("maintain", me[which(names(me)==1)], me[which(names(me)==0)]), sep="__")
						}
					}
					else{
						if(which(names(me)==1)> which(names(me)==0)){
							out[tt]="up"
						}
						if(which(names(me)==1)< which(names(me)==0)){
							out[tt]="down"
						}
						if(which(names(me)==1)== which(names(me)==0)){
							out[tt]="maintain"
						}
					}
					names(out)=sel_pair
					
					return(out)
				}
			}
			
			out=sapply(1:length(com_sel_by), fun)
			out=t(out)
			#rownames(out)=genename
			out=cbind(genename, category, out)
			return(out)
		}
	}
	
	####omics data separate
	mrna=omics[[which(category=="mrna")]]
	omics_new<<-omics[-which(category=="mrna")]
	category_new<<-category[-which(category=="mrna")]
	
	
	####pathway level
	es.max <<- gsva(mrna, pathway_list, mx.diff=FALSE, verbose=FALSE, parallel.sz=cores, method=gsva_method)
	rownames(es.max)=names(pathway_list)

	out<<-pathway_list
	name<<-names(pathway_list)

	core_all=function(path_name){
		##Step_01: find the alter matrix
		shan=path_alter_find(path_name)
		#save(shan, file="shan.Rdata")
		score=es.max[which(rownames(es.max)==path_name),]
		score=cbind(names(score), score)

		##Step_02: compare the alter
		com=score_compar(alter=shan, score=score, type=category_new, split_by="!!@@")

		##Step 03: merge the outcome of the compare
		return(path_driver_find(dat=com, q_cutoff=0.1, sel_pair=c("1--0", "-1--0", "-1--1"), ps_cutoff=0.01, split_by="--", return_mean=TRUE))
	}

	cl <<- makeCluster(cores); clusterExport(cl, c("out", "name", "omics_new", "category_new", "es.max", "path_alter_find", "score_compar", "path_driver_find"))#; clusterEvalQ(cl, "DriverFind")
	see=parLapply(cl, name, core_all)
	stopCluster(cl)

	names(see)=name

	filter_see=function(x){
		if (is.null(dim(x))){
			if (all(is.na(x))){
				return(x)
			}
		}
		else{
			tag=intersect(which(as.character(as.matrix(x[,3]))=="-"), which(as.character(as.matrix(x[,4]))=="-"))
			if (length(tag)==0){
				return(x)
			}
			else{
				return(x[-tag, ])
			}
		}
	}

	see1=lapply(see, filter_see)
	
	return(see1)
}



analysis=function(x){
	tit=names(x)
	x=x[[1]]
	
	prepare=function(i){
		you=data.table::as.data.table(x[,i]); colnames(you)="V1"
		you=splitstackshape::cSplit(you, "V1", "__")
		
		if (ncol(you)!=3){
			return(NULL)
		}
		if (ncol(you==3)){
			###calculate mean value
			me=na.omit(you); me=as.numeric(as.matrix(me[,2]))-as.numeric(as.matrix(me[,3]))
			
			###plot data
			you=cbind(x[which(as.character(as.matrix(you[,1]))!="-"),1], x[which(as.character(as.matrix(you[,1]))!="-"),2], me)
			you=as.data.frame(you); colnames(you)=c("V1", "V2", "V3")
			you[,1]=as.character(as.matrix(you[,1])); you[,2]=as.character(as.matrix(you[,2])); you[,3]=as.numeric(as.matrix(you[,3]))
			you=you[order(you[,3]), ]
			return(you)
		}
	}
	
	cna_col="#C51A7D"; cnd_col="#B8E085"; meth_col="#92C5DE"; mut_col="yellow"
	
	if (nrow(x)==0){
		stop("No potential driver genomic alterations has been identified.")
	}
	if (is.null(dim(x))){
		x=t(as.data.frame(x))
		tag=as.character(as.matrix(x[,3:4]))
		sel=which(tag!="-")
		if (length(sel)==1){
			me=data.table::as.data.table(tag[sel]); colnames(me)="V1"; me=splitstackshape::cSplit(me, "V1", "__")
			me=as.numeric(as.matrix(me[,2]))-as.numeric(as.matrix(me[,3]))
			if (sel==1){
				if (as.character(as.matrix(x[,2]))=="cnv"){
					pl=t(as.data.frame(c(as.character(as.matrix(x[,1])), "cna", me)))
				}
				else{
					pl=t(as.data.frame(c(as.character(as.matrix(x[,1])), as.character(as.matrix(x[,2])), me)))
				}
			}
			else{
				pl=t(as.data.frame(c(as.character(as.matrix(x[,1])), "cnd", me)))
			}
		}
		else{
			for (i in 1:length(sel)){
				me=data.table::as.data.table(tag[sel[i]]); colnames(me)="V1"; me=splitstackshape::cSplit(me, "V1", "__")
				me=as.numeric(as.matrix(me[,2]))-as.numeric(as.matrix(me[,3]))
				if (i==1){
					if (as.character(as.matrix(x[,2]))=="cnv"){
						pl=t(as.data.frame(c(as.character(as.matrix(x[,1])), "cna", me)))
					}
					else{
						pl=t(as.data.frame(c(as.character(as.matrix(x[,1])), as.character(as.matrix(x[,2])), me)))
					}
				}
				else{
					pl=rbind(pl, t(as.data.frame(c(as.character(as.matrix(x[,1])), "cnd", me))))
				}
			}
		}
	}
	else{
		if (nrow(x)==1){
			sels=which(as.character(as.matrix(x[3:4]))!="-")
			if (length(sels)==1){
				me=data.table::as.data.table(x[,c(3:4)[sels]]); colnames(me)="V1"; me=splitstackshape::cSplit(me, "V1", "__")
				me=as.numeric(as.matrix(me[,2]))-as.numeric(as.matrix(me[,3]))
				pl=c(as.character(as.matrix(x[,1:2])), me)
				pl=t(as.data.frame(pl)); colnames(pl)=c("V1", "V2", "V3")
				
				typp=as.character(as.matrix(pl[,2]))
				if (typp=="cnv"){
					if (sels==1){
						pl[,2]="cna"
					}
					else{
						pl[,2]="cnd"
					}
				}
			}
			else{
				pl=alist()
				for (i in sels){
					me=data.table::as.data.table(x[,c(3:4)[sels]]); colnames(me)="V1"; me=splitstackshape::cSplit(me, "V1", "__")
					me=as.numeric(as.matrix(me[,2]))-as.numeric(as.matrix(me[,3]))
					pl[[i]]=t(as.data.frame(c(as.character(as.matrix(x[,1:2])), me)))
					typp=as.character(as.matrix(pl[[i]][,2]))
					if (typp=="cnv"){
						if (sels==1){
							pl[[i]][,2]="cna"
						}
						else{
							pl[[i]][,2]="cnd"
						}
					}
				}
				for (i in 1:length(pl)){
					if (i==1){
						sl=pl[[1]]
					}
					else{
						sl=rbind(sl, pl[[1]])
					}
				}
				pl=sl; pl=as.data.frame(pl); colnames(pl)=c("V1", "V2", "V3")
			}
		}
		if (nrow(x)>1){
			you=prepare(3)
			zuo=prepare(4)
			if (!is.null(zuo)){
				zuo[which(zuo[,2]=="cnv"),2]="cnd"
			}
			if (!is.null(you)){
				you[which(you[,2]=="cnv"),2]="cna"
			}
			pl=as.data.frame(rbind(zuo, you))
		}
	}
	cate=names(table(pl[,2]))
	for (i in 1:length(cate)){
		if (i==1){
			pl1=pl[which(pl[,2]==cate[i]),]
		}
		else{
			pl1=rbind(pl1, pl[which(pl[,2]==cate[i]),])
		}
	}
		
	pl=pl1
	if (all(is.null(dim(pl)), length(pl)>0)){
		pl=t(as.data.frame(pl))
	}
	pl[which(duplicated(pl[,1])),1]=paste0(pl[which(duplicated(pl[,1])),1], "_", pl[which(duplicated(pl[,1])),2])
	library(ggplot2)
	pl=as.data.frame(pl); colnames(pl)=c("V1", "V2", "V3")
	pl$V1 <- factor(pl$V1,levels=as.character(as.matrix(pl$V1)), ordered=TRUE)
	
	all_col=c(cna_col, cnd_col, meth_col, mut_col); names(all_col)=c("cna", "cnd", "met", "mut")
		
	cc=match(names(table(pl$V2)), names(all_col))
	for (i in 1:length(cc)){
		if (i==1){
			col=rep(all_col[cc[i]], table(pl$V2)[i])
		}
		else{
			col=c(col, rep(all_col[cc[i]], table(pl$V2)[i]))
		}
	}
		
	p=ggplot(data = pl, mapping = aes(x = V1, y = V3, fill=V2)) + 
	scale_fill_manual(values = col)+ 
	geom_bar(stat= 'identity') + theme(axis.text.x = element_text(angle = 45, hjust =  1))+
	xlab("genes")+ ylab("difference of mean value")+ guides(fill=guide_legend(title=NULL))+
	labs(title = tit)+ theme(plot.title= element_text(face= "bold", color="black", vjust=0.5, hjust=0.5, angle=0))
	
	re=alist()
	re[[1]]=p; re[[2]]=pl
	return(re)
}
