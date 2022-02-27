# Lagrangian for planar orbit in sperical coordinates

For a bodies of mass M and m, where M >> m, and r is the distance of m from M,
the lagrangian can be written as

$$
\mathscr{L} = \frac{1}{2}(m \dot{r}^2 + m r^2 \dot{\theta}^2) + \frac{GMm}{r}
$$

Some useful partial derivatives to substitute in Euler-Lagrange equations

$$
\frac{\partial \mathscr{L}}{\partial r} = m r \dot{\theta}^2 - \frac{GMm}{r^2}
$$

$$
\frac{\partial \mathscr{L}}{\partial \dot{r}} = m \dot{r}
$$

$$
\frac{\partial \mathscr{L}}{\partial \theta} = 0
$$

$$
\frac{\partial \mathscr{L}}{\partial \dot{\theta}} = m r^2 \dot{\theta}
$$

Simplifying Euler-Lagrange equation for the r coordinate

$$
0 = \frac{\partial \mathscr{L}}{\partial r} -
        \frac{d}{dt}(\frac{\partial \mathscr{L}}{\partial \dot{r}}) = r
        \dot{\theta}^2 - \frac{GM}{r^2} -\ddot{r}
$$

$$
\longrightarrow \ddot{r} = r \dot{\theta}^2 - \frac{GM}{r^2}
$$

Simplifying Euler-Lagrange equation for the $$\theta$$ coordinate

$$
0 = \frac{\partial \mathscr{L}}{\partial \theta} -
        \frac{d}{dt}(\frac{\partial \mathscr{L}}{\partial \dot{\theta}}) =
        \frac{d}{dt}(m r^2 \dot{\theta})
$$

Simplifying and integrating we get

$$
\longrightarrow L = m r^2 \dot{\theta}
$$

where L is the angular momentum

Substituting for $$\dot{\theta}$$ we get

$$
\ddot{r} = \frac{L^2}{m^2r^3} - \frac{GM}{r^2}
$$
