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

positive variables

  plant(crop)
  purchase(need, scenario)
  sell(surplus, scenario)
        
variables

*** model variables

  z                      current objective value

*** calculation variables

  WSobjectives(scenarios)  objective values for each WS scenario
  EVobjectives(scenarios)  objective values for each EV scenario
  WS                       wait-and-see value
  EVV                      expected result of using the expected value solution

equations

  obj                         minimize negative profit
  land                        must not exceed maximum land area
  quotaWheat(scenario)        wheat quota must be met
  quotaCorn(scenario)         corn quota must be met
  quotaBeets(scenario)        beets must not exceed quota
  quotaExcessBeets(scenario)  excessBeets exceed quota;
  
obj..
  z
  =E= sum(crop, plantingCost(crop) * plant(crop))
      - sum(scenario, scenarioProb(scenario)
        * (sum(surplus, sell(surplus, scenario) * salePrice(surplus))
           - sum(need, purchase(need, scenario) * purchasePrice(need))));

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
  
model hw1_5a /all/;
option solprint = off;

*** compute the wait-and-see objective value WS

loop(scenarios,
  loop(crop,
    yield("1", crop) = yields(scenarios, crop);
    )
  solve hw1_5a using lp minimizing z;
  WSobjectives.l(scenarios) = z.l;
  );
WS.l =  sum(scenarios, scenariosProb(scenarios)*WSobjectives.l(scenarios));

*** compute the expected value of the yields

loop(crop,
  yield("1", crop) = sum(scenarios, scenariosProb(scenarios)*yields(scenarios, crop));
  );

*** find the expected value solution

solve hw1_5a using lp minimizing z;

*** compute the expected result of using the expected value solution EEV

loop(crop,
  plant.fx(crop) = plant.l(crop);
  )

loop(scenarios,
  loop(crop,
    yield("1", crop) = yields(scenarios, crop);
    )
  solve hw1_5a using lp minimizing z;
  display plant.l;
  EVobjectives.l(scenarios) = z.l;
  );
EVV.l =  sum(scenarios, scenariosProb(scenarios)*EVobjectives.l(scenarios));

*** output data

file out / farmersPH.put /
put out;
put 'Wait-and-see objective value';
put WS.l /;
put 'Expected result of using the expected value solution';
put EVV.l;
putclose out;
