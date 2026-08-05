Your ROCm version is older then 7.0.2.
From our measurements, earlier ROCm versions offer
underwhelming performance on eigensolvers, that you will
notice if you run Ground-State calculations.

Performance hit is especially bad when using a ROCm older
than 6.4.2, hence we STRONGLY advise you update to a more
recent ROCm version in that case.

ROCm versions between 6.4.2 and 7.0.2 show acceptable
performance but may consume too much memory on big cases.
