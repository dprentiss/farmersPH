* David Prentiss
* Farmer's problem with progressive hedging

sets

  crop       crops we can plant        / wheat, corn, beets /
  need       crops we may need to buy  / wheat, corn /
  surplus    crops we can sell         / wheat, corn, beets, excessBeets / 
  scenario   current yield scenario    / 1 /
  scenarios  yields for loop           / 1*9 /

parameters

  plantingCost(crop)        / wheat 150
                              corn  230
                              beets 260 /

  purchasePrice(need)       / wheat 283
                              corn  210 /

  salePrice(surplus)        / wheat       170
                              corn        150
                              beets       36
                              excessBeets 10 /

  minimum(need)             / wheat 200
                              corn  240 /

  maxLand                   / 500 /
  
  beetQuota                 / 6000 /

*** to reuse the RP code we give a probability of
*   one to each wait-and-see scenario

  scenarioProb(scenario)    / 1 1 /
  
*** probabilities for each scenario

  scenariosProb(scenarios)  / 1 0.2222222222222222
                              2 0.1111111111111111
                              3 0.1111111111111111
                              4 0.0555555555555555
                              5 0.0555555555555555
                              6 0.0555555555555555
                              7 0.0555555555555555
                              8 0.1111111111111111
                              9 0.2222222222222222 /

*** all yields by scenario

table yields(scenarios, crop)

   wheat   corn    beets

1  2.0     2.4     16
2  2.0     2.4     20 
3  2.0     2.4     24
4  2.5     3.0     16
5  2.5     3.0     20
6  2.5     3.0     24
7  3.0     3.6     16
8  3.0     3.6     20
9  3.0     3.6     24
   ;

*** current scenario yields

table yield(scenario, crop)

   wheat   corn    beets

1  0       0       0 
   ;

scalar rho / 3 /;
scalar threshold / 0.1 /

positive variables

  plant(crop)
  purchase(need, scenario)
  sell(surplus, scenario)
        
variables

*** model variables

  z                      current objective value

*** calculation variables

  ExpPlant(crop)
  g
  ws(scenarios, crop)
  w(scenario, crop)

equations

  obj0                        minimize negative profit
  objk                        minimize negative profit
  land                        must not exceed maximum land area
  quotaWheat(scenario)        wheat quota must be met
  quotaCorn(scenario)         corn quota must be met
  quotaBeets(scenario)        beets must not exceed quota
  quotaExcessBeets(scenario)  excessBeets exceed quota;
  
obj0..
  z
  =E= sum(crop, plantingCost(crop) * plant(crop))
      - sum(scenario, scenarioProb(scenario)
        * (sum(surplus, sell(surplus, scenario) * salePrice(surplus))
           - sum(need, purchase(need, scenario) * purchasePrice(need))));
        
objk..
  z
  =E= sum(crop, plantingCost(crop) * plant(crop))
      - sum(scenario, scenarioProb(scenario)
        * (sum(surplus, sell(surplus, scenario) * salePrice(surplus))
           - sum(need, purchase(need, scenario) * purchasePrice(need))))
      + sum((scenario, crop), w(scenario, crop) * plant(crop))
      + rho/2 * sqrt(sum(crop, (plant.l(crop) - ExpPlant.l(crop))**2))**2;

land..
  sum(crop, plant(crop))
  =L= maxLand;
  
quotaWheat(scenario)..
  yield(scenario, "wheat" ) * plant("wheat") 
    + purchase("wheat", scenario)
    - sell("wheat", scenario)
  =G= minimum("wheat");
  
quotaCorn(scenario)..
  yield(scenario, "corn") * plant("corn") 
    + purchase("corn", scenario)
    - sell("corn", scenario)
  =G= minimum("corn");

quotaBeets(scenario)..
  sell("beets", scenario)
  =L= beetQuota;

quotaExcessBeets(scenario)..
  sell("beets", scenario) + sell("excessBeets", scenario)
  =L= yield(scenario, "beets") * plant("beets");
  
model firstIter / obj0, land, quotaWheat, quotaCorn, quotaBeets, quotaExcessBeets /;
model kthIter / objk, land, quotaWheat, quotaCorn, quotaBeets, quotaExcessBeets /;
option solprint = off;

*** init output

file out / farmersPH.put /
put out;

*** solve each scenario without penalties

loop(crop,
  ExpPlant.l(crop) = 0;
  );
loop(scenarios,
  loop(crop,
    yield("1", crop) = yields(scenarios, crop);
    );
  solve firstIter using lp minimizing z;
  loop(crop,
    ExpPlant.l(crop) = ExpPlant.l(crop) + scenariosProb(scenarios) * plant.l(crop);
    ws.l(scenarios, crop) = rho * (plant.l(crop) - ExpPlant.l(crop));
    );
  );

put '1'; loop(crop, put ExpPlant.l(crop)); put /;

g.l = 1;
while(g.l gt threshold,
  loop(crop,
    ExpPlant.l(crop) = 0;
    );
  g.l = 0;
  loop(scenarios,
    loop(crop,
      w.l("1", crop) = ws.l(scenarios, crop);
      yield("1", crop) = yields(scenarios, crop);
      );
    solve kthIter using nlp minimizing z;
    loop(crop,
      ExpPlant.l(crop) = ExpPlant.l(crop) + scenariosProb(scenarios) * plant.l(crop);
      ws.l(scenarios, crop) = ws.l(scenarios,  crop) + rho * (plant.l(crop) - ExpPlant.l(crop));
      );
    );
  put 'k'; loop(crop, put ExpPlant.l(crop)); put /;
  );

*** close output
putclose out;
