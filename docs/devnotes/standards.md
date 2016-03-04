The following practices should be followed in gPhoton unless there is a well articulated reason to do otherwise (that is explained by notes in code):

* All occurrences of `UNION` in SQL should be `UNION ALL.`
* All ranges should be inclusive of the lower boundary and exclusive of the upper boundary. e.g. `t0 <= t < t1`
