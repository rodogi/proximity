# Proximity

The proximity distance $d$ between two sets of nodes $S$ and $T$ is defined as

$$
d(S, T) = \frac{1}{||T||}\sum_{t \in T} \min_{s \in S} \delta(s, t),
$$

where $\delta(s, t)$ is the shortest path length between nodes $s$ and $t$.