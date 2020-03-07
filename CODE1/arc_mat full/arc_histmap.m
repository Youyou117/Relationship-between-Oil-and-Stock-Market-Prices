function arc_histmap(pltvariables,results,options)
    data = am_DataSource(pltvariables,results,options);
    hist = am_Histogram(data);
    map = am_BaseMap(data);
    link = am_Linkage;
    link.addMap(map);
    link.addGraph(hist);
end
