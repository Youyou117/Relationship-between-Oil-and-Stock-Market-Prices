function arc_sarmap(y,x,latt,long,results,options)
    data = am_DataSource(y,results,options);
    lsar = am_LocalSAR(data,y,x);
    gsar = am_GlobalSAR(data,y, x);
    map = am_BaseMap(data);
    link = am_Linkage;
    link.addMap(map);
    link.addGraph(lsar);
    link.addGraph(gsar);
end