* David Prentiss
* Farmer's problem with progressive hedging

sets

  crop       crops we can plant        / wheat, corn, beets /
  need       crops we may need to buy  / wheat, corn /
  surplus    crops we can sell         / wheat, corn, beets, excessBeets /
  scenarios  yields for loop           / 1*3 /

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

  yield(crop)
  prob
  Eplant(crop)
  g
  ws(scenarios, crop)
  w(crop)

*** probabilities for each scenario
  scenarioProb(scenarios)   / 1 0.3333333333333333
                              2 0.3333333333333333
                              3 0.3333333333333333 /

*** all yields by scenario
table yields(scenarios, crop)
   wheat   corn    beets
1  2.0     2.4     16
2  2.5     3.0     20
3  3.0     3.6     24

scalar
  rho / 2 /
  threshold / 0.1 /

positive variables
  plant(crop)
  purchase(need)
  sell(surplus)

variables
  z                           current objective value

equations

  obj0                        minimize negative profit
  objk                        minimize negative profit
  land                        must not exceed maximum land area
  quotaWheat                  wheat quota must be met
  quotaCorn                   corn quota must be met
  quotaBeets                  beets must not exceed quota
  quotaExcessBeets            excessBeets exceed quota
  ;

obj0..
  z
  =E=
  sum(crop, plantingCost(crop) * plant(crop))
      - prob * (sum(surplus, sell(surplus) * salePrice(surplus))
           - sum(need, purchase(need) * purchasePrice(need)));

objk..
  z
  =E=
  sum(crop, plantingCost(crop) * plant(crop))
      - prob * (sum(surplus, sell(surplus) * salePrice(surplus))
           - sum(need, purchase(need) * purchasePrice(need)))
      + sum(crop, w(crop) * plant(crop))
      + rho/2 * sqrt(sum(crop, power((plant(crop) - EPlant(crop)), 2)))**2;

land..
  sum(crop, plant(crop))
  =L= maxLand;

quotaWheat..
  yield("wheat") * plant("wheat") + purchase("wheat") - sell("wheat")
  =G= minimum("wheat");

quotaCorn..
  yield("corn") * plant("corn") + purchase("corn") - sell("corn")
  =G= minimum("corn");

quotaBeets..
  sell("beets")
  =L= beetQuota;

quotaExcessBeets..
  sell("beets") + sell("excessBeets")
  =L= yield("beets") * plant("beets");

model firstIter / obj0, land, quotaWheat, quotaCorn, quotaBeets, quotaExcessBeets /;
model kthIter / objk, land, quotaWheat, quotaCorn, quotaBeets, quotaExcessBeets /;
option solprint = off;

*** init output
file out / farmersPH.put /;
put out;

*** solve each scenario without penalties
loop(crop, Eplant(crop) = 0;);
loop(scenarios,
  prob = scenarioProb(scenarios);
  loop(crop, yield(crop) = yields(scenarios, crop););
  solve firstIter using lp minimizing z;
  loop(crop,
    Eplant(crop) = Eplant(crop) + prob * plant.l(crop);
    ws(scenarios, crop) = rho * (plant.l(crop) - EPlant(crop));
    );
  );

put '1'; loop(crop, put EPlant(crop)); put /;

*** solve penalized version until converged
g = .5;
while(g gt threshold,
  loop(crop, EPlant(crop) = 0;);
  loop(scenarios,
    prob = scenarioProb(scenarios);
    loop(crop,
      w(crop) = ws(scenarios, crop);
      yield(crop) = yields(scenarios, crop);
      );
    solve kthIter using nlp minimizing z;
    loop(crop,
      EPlant(crop) = EPlant(crop) + prob * plant.l(crop);
      ws(scenarios, crop) = ws(scenarios, crop) + rho * (plant.l(crop) - EPlant(crop));
      );
    );
    g = g - 0.1;
  put 'k'; loop(crop, put EPlant(crop)); put /;
  );

*** close output
putclose out;
