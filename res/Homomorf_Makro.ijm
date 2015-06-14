T = getDirectory("Choose a Directory")
T = T+"properties";
run("RiceHommomorfEst ", "choose=&T");
selectWindow("rician");
selectWindow("gauss");
run("Jet");
selectWindow("rician");
run("Jet");
