h=0.5
dO=-0.5

for(dY in seq(0.05,0.45,0.1)) {
for(so in seq(0.1,0.9,0.1)) {
for(ssro in c(0,0.2,0.5)) {
  if(ssro>1) ssro=1  
for(ssry in seq(0,0.5,0.1)) {
for(ssr in seq(0.1,0.9,0.05)) {
for(u in c(0,0.001)) {
  for(lambda in 10^(seq(0,3,3))) {
  
  xf=0.999
  xm=0.999
  y=1-u/so

  t=0

  xflist=c(xf)
  xmlist=c(xm)
  ylist=c(y)
  
  while(t<20000) {
    
    wbarxf=xf*xm+(1-h*ssr)*(xm*(1-xf)+xf*(1-xm))+(1-ssr)*(1-xf)*(1-xm);
    wbarxm=1/2*xf*y+1/2*(1-so)*xf*(1-y)+(1/2+dY)*y*(1-xf)*(1-ssry)+(1/2+dO)*(1-ssro)*(1-y)*(1-xf);
    wbary=1/2*xf*y+1/2*(1-so)*xf*(1-y)+(1/2-dY)*y*(1-xf)*(1-ssry)+(1/2-dO)*(1-ssro)*(1-y)*(1-xf);
    
    xfnext=(xf*xm+1/2*(1-h*ssr)*(xm*(1-xf)+xf*(1-xm)))/wbarxf
    xmnext=(1/2*xf*y+1/2*(1-so)*xf*(1-y))/wbarxm
    ynext=(1/2*xf*y*(1-u)+(1/2-dY)*(1-xf)*y*(1-u*lambda)*(1-ssry))/wbary
    
    xf=xfnext
    xm=xmnext
    y=ynext
    
    xflist=append(xflist,xf)
    xmlist=append(xmlist,xm)
    ylist=append(ylist,y)
    
    t=t+1
  }
  
  
  xfsrlist=1-xflist
  olist=1-ylist
  xmsrlist=1-xmlist
  
  cat(t,u,lambda,u*lambda,dY,dO,ssr,ssry,ssro,so,h,1-xf,1-xm,1-y,u/so,
      if(var(tail(olist,2000))<0.0001) 1 else 0,
      if(var(tail(xmsrlist,2000))<0.0001) 1 else 0,
      if(var(tail(xfsrlist,2000))<0.0001) 1 else 0,
      if(max(olist)<=0.01) max(olist) else 1,
      if(max(xmsrlist)<=0.01) max(xmsrlist) else 1,
      if(max(xfsrlist)<=0.01) max(xfsrlist) else 1,sep=",")
  cat("\n")
  
}
}
}
}
}
}
}