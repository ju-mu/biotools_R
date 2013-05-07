
#' Calculates donor side shapiro scores 
#' @param ref_fasta DNAStringSet or list of strings with the reference sequences
#' @param seq_fasta DNAStringSet or list of strings with the target sequences
#' @returnType numeric vector
#' @return vector with Shapiro scores
#' 
#' @author Julius Muller
#' @export
shapiro_donor<-function(ref_fasta,seq_fasta){
  rdict<-list("A"=1,"C"=2,"G"=3,"T"=4)
  
  if(class(ref_fasta)!="DNAStringSet")ref_fasta<-readDNAStringSet(ref_fasta)
  
  mpwm<-consensusMatrix(ref_fasta,as.prob=T,baseOnly=T)[1:4,]*100
  
  tmin<-apply(mpwm,2,which.min)
  tmin<-sum(mpwm[cbind(tmin,seq(tmin))])
  tmax<-apply(mpwm,2,which.max)
  tmax<-sum(mpwm[cbind(tmax,seq(tmax))])
  
  if(class(seq_fasta)=="DNAStringSet"){dss2<-as.character(seq_fasta)
  }else{dss2<-as.character(readDNAStringSet(seq_fasta))}
  dss2<-lapply(dss2,strsplit,"")
  
  unlist(lapply(dss2,function(x){
    mpos<-unlist(lapply(unlist(x),function(y){rdict[[y]]}))
    100*(sum(mpwm[cbind(mpos,seq(mpos))])-tmin)/(tmax-tmin)
  }))
}

#' Calculates acceptor side shapiro scores 
#' @param ref_fasta DNAStringSet or list of strings with the reference sequences
#' @param seq_fasta DNAStringSet or list of strings with the target sequences
#' @returnType numeric vector
#' @return vector with Shapiro scores
#' 
#' @author Julius Muller
#' @export
shapiro_acceptor<-function(ref_fasta,seq_fasta){
  rdict<-list("A"=1,"C"=2,"G"=3,"T"=4)
  
  names(seq_fasta)<-1:length(seq_fasta)
  
  if(class(ref_fasta)!="DNAStringSet")ref_fasta<-readDNAStringSet(ref_fasta)
  if(class(seq_fasta)=="DNAStringSet"){ospyri<-as.character(seq_fasta)
  }else{ospyri<-as.character(readDNAStringSet(seq_fasta))}
  
  mpwm<-consensusMatrix(ref_fasta,as.prob=T,baseOnly=T)[1:4,]*100
  stopifnot(all(dim(mpwm)==c(4,15)))
  
  pyri<-mpwm[,1:10]
  cagg<-mpwm[,11:14]
  
  l1<-apply(pyri,2,which.min)
  l1<-pyri[cbind(l1,seq(l1))]
  l1<-sum(sort(l1,decreasing=F)[1:8])
  l2<-apply(cagg,2,which.min)
  l2<-sum(cagg[cbind(l2,seq(l2))])
  
  h1<-apply(pyri,2,which.max)
  h1<-pyri[cbind(h1,seq(h1))]
  h1<-sum(sort(h1,decreasing=T)[1:8])
  h2<-apply(cagg,2,which.max)
  h2<-sum(cagg[cbind(h2,seq(h2))])
  
  ospyri<-lapply(ospyri,strsplit,"")
  spyri<-lapply(ospyri,function(x){x[[1]][1:10]})
  t1<-lapply(spyri,function(x){
    mpos<-unlist(lapply(unlist(x),function(y){rdict[[y]]}))
    sum(sort(pyri[cbind(mpos,seq(mpos))],decreasing=T)[1:8])
  })
  scagg<-lapply(ospyri,function(x){x[[1]][11:14]})
  t2<-lapply(scagg,function(x){
    mpos<-unlist(lapply(unlist(x),function(y){rdict[[y]]}))
    sum(cagg[cbind(mpos,seq(mpos))])
  })
  
  fscores<-lapply(names(t1),function(x){
    100*((t1[[x]]-l1)/(h1-l1)+(t2[[x]]-l2)/(h2-l2))/2
  })
  names(fscores)<-names(t1)
  unlist(fscores)
}
