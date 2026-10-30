
vecmseinput(simul) = norm(w-Kv)^2/n;
vecmseitdepHTminmse(simul) = norm(v1itdepHTminmse-v)^2/m;
vecmseitdepHTmingcv(simul) = norm(v1itdepHTmingcv-v)^2/m;
vecmsefixedHT(simul) = norm(v1fixedHT-v)^2/m;
vecmseitdepSTminmse(simul) = norm(v1itdepSTminmse-v)^2/m;
vecmseitdepSTmingcv(simul) = norm(v1itdepSTmingcv-v)^2/m;
vecmsefixedST(simul) = norm(v1fixedST-v)^2/m;

vecfalsepositdepHTminmse(simul) = ...
                              sum((abs(v1itdepHTminmse)>eps).*(abs(v)<eps));
vecfalsepositdepHTmingcv(simul) = ...
                              sum((abs(v1itdepHTmingcv)>eps).*(abs(v)<eps));
vecfalseposfixedHT(simul) = sum((abs(v1fixedHT)>eps).*(abs(v)<eps));
vecfalsepositdepSTminmse(simul) = ...
                              sum((abs(v1itdepSTminmse)>eps).*(abs(v)<eps));
vecfalsepositdepSTmingcv(simul) = ...
                              sum((abs(v1itdepSTmingcv)>eps).*(abs(v)<eps));
vecfalseposfixedST(simul) = sum((abs(v1fixedST)>eps).*(abs(v)<eps));

vecfalsenegitdepHTminmse(simul) = ...
                              sum((abs(v1itdepHTminmse)<eps).*(abs(v)>eps));
vecfalsenegitdepHTmingcv(simul) = ...
                              sum((abs(v1itdepHTmingcv)<eps).*(abs(v)>eps));
vecfalsenegfixedHT(simul) = sum((abs(v1fixedHT)<eps).*(abs(v)>eps));
vecfalsenegitdepSTminmse(simul) = ...
                              sum((abs(v1itdepSTminmse)<eps).*(abs(v)>eps));
vecfalsenegitdepSTmingcv(simul) = ...
                              sum((abs(v1itdepSTmingcv)<eps).*(abs(v)>eps));
vecfalsenegfixedST(simul) = sum((abs(v1fixedST)<eps).*(abs(v)>eps));

vectruepositdepHTminmse(simul) = ...
                              sum((abs(v1itdepHTminmse)>eps).*(abs(v)>eps));
vectruepositdepHTmingcv(simul) = ...
                              sum((abs(v1itdepHTmingcv)>eps).*(abs(v)>eps));
vectrueposfixedHT(simul) = sum((abs(v1fixedHT)>eps).*(abs(v)>eps));
vectruepositdepSTminmse(simul) = ...
                              sum((abs(v1itdepSTminmse)>eps).*(abs(v)>eps));
vectruepositdepSTmingcv(simul) = ...
                              sum((abs(v1itdepSTmingcv)>eps).*(abs(v)>eps));
vectrueposfixedST(simul) = sum((abs(v1fixedST)>eps).*(abs(v)>eps));
