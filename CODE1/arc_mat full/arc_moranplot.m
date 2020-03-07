function arc_moranplot(variable,W,results,options)
    data = am_DataSource(variable,results,options);
    moran = am_Moranplot(data);
    map = am_BaseMap(data);
    link = am_Linkage;
    link.addMap(map);
    link.addGraph(moran);
end
