This is a post about rotation operations, that is operations that modify an orientation (or attitude). There are a few ways we can represent rotations: unit quaternions, rotation matrices and axis-angle representations are some of the most useful. Rotation matrices and unit quaternions are also commonly used to represent the orientations themselves, and axis-angle is most common when talking about changes in orientation (it's basically an "angular velocity" representation).

There are two things we might want to do with a given rotation operator: rotate a single vector, or rotate an entire orientation. When we're in quaternion representation, rotating an orientation looks like simple (quaternion) multiplication:
$$
q_1 = \Delta{q} \ q_0
$$
Similarly, for rotation matrices we can use matrix multiplication:
$$
R_1 = \Delta{R}\ R_0
$$
When rotating a vector, rotation matrices still work by matrix multiplication (or matrix-vector multiplication, if you consider that different):
$$
v_1 = \Delta{R}\ v_0
$$
However, quaternions use an operation known as "conjugation":
$$
v_1 = \Delta{q}\ v_0\ \Delta{q}^{-1} = \Delta{q}\ v_0\ \Delta{q}^*,
$$
where the second equality follows from the equality of the quaternion inverse (reciprocal) and the quaternion conjugate (computed by negating its vector part) for unit quaternions. The first holds even for non-unit quaternions, though I will not use that fact here.

Typically, the way that these operations are performed for axis-angle rotations is to first compute a corresponding rotation matrix or quaternion, then apply one of the above procedures. This can be an efficient approach if the derived matrix or quaternion is going to be used many times, but doing so once and then throwing it away can be wasteful.

First, I'll look at rotating a vector according to an angular velocity vector $\omega$ and a delta-time scale factor $\Delta{t}$. The total scaled-axis representation of the rotation is then $\omega\Delta{t}$, with total angle $\phi = |\omega| \Delta{t}$ and axis $\hat{n} = \omega / |\omega|$ the unit vector parallel to the angular velocity. This is a common representation in applications. 

The quaternion form of this scaled-axis rotation has real (scalar) part $\cos(\omega\Delta{t}/2)$ and imaginary (vector) part $\hat{n} \sin(\omega \Delta{t}/2)$. Because the variables we have to work with are the angular velocity and the time step, *not* the unit axis $\hat{n}$, computing this requires normalising the angular velocity (involving a square root) as well as computing a sine and cosine. Not ideal.

How about the rotation matrix representation? Rodrigues' rotation formula gives us:
$$
R(\omega\Delta{t}) = \mathbb{I} + \sin(|\omega|\Delta{t}) [\hat{n}]_\times + (1 - \cos(|\omega|\Delta{t})) [\hat{n}]_\times^2
$$
Here, $\mathbb{I}$ is the identity matrix and the antisymmetric matrix $[a]_\times$ is such that $[a]_\times \cdot v = a \times v$ is the usual cross product. Therefore, the effect on a vector $v$ is:
$$
R(\omega\Delta{t}) \cdot v = v + \sin(|\omega|\Delta{t}) (\hat{n} \times v) + (1 - \cos(|\omega|\Delta{t})) \hat{n} \times (\hat{n} \times v)
$$
This also involves normalising the angular velocity and computing some trig functions. Can we do better?

One approach is "first order and normalise": take the rotation matrix to first order in $\Delta{t}$ and then normalise the resulting vector to the same length as the input:
$$
\begin{align}
v' &= v_0 + (|\omega| \Delta{t}) (\hat{n} \times v_0) = v_0 + (\omega \times v_0) \Delta{t} \\
v_1 &= \sqrt\frac{|v_0|^2}{|v'|^2} v' = \frac{v_0 + (\omega \times v_0)\Delta{t}}{\sqrt{1 + (|\omega\times v_0|^2 / v_0^2) \Delta{t}^2}}
\end{align}
$$
This looks okay at first glance - now just a square root is necessary and no trig functions, but there's a problem: This doesn't actually rotate around $\hat{n}$. Notice that $\omega$ only appears in the combination $\omega \times v_0$, which is the same as $(\omega + \alpha v_0)\times v_0$ for any $\alpha$, because cross products between a vector and itself are zero. This operation actually rotates around $\omega - (\omega \cdot v_0 / v_0^2) v_0$, orthogonal to $v_0$. 

It's notable that this same sort of problem doesn't occur in the quaternion domain: The first-order approximation to the quaternion is:
$$
q'(\omega \Delta{t}) = (\cos(|\omega|\Delta{t}/2), \hat{n}\sin(|\omega|\Delta{t}/2)) \approx (1, \hat{n}(|\omega|\Delta{t}/2)) = (1, \omega \Delta{t}/2).
$$
Dividing by the magnitude of this denormalised quaternion gives:
$$
q_1 = \frac{[1, \omega\Delta{t}/2]}{\sqrt{1 + \omega^2 \Delta{t}^2 / 4}}
$$
which has its imaginary vector component in the same direction as $\omega$. It's not too hard to figure out the angle this quaternion turns through:
$$
\begin{align}
\cos(\Theta/2) &= \frac{1}{\sqrt{1 + \omega^2 \Delta{t}^2 / 4}} \\
\cos(\Theta/2)^{-2} &= 1 + \omega^2 \Delta{t}^2 /4 \\
\left(\frac{\omega \Delta{t}}{2}\right)^2 &= \sec^2(\Theta/2) - 1 = \tan^2(\Theta/2) \\
\Theta &= 2 \arctan(|\omega|\Delta{t}/2) \approx |\omega| \Delta{t}
\end{align}
$$
Basically, the angle saturates according to the arctangent of the desired angle, leading to absolute angles between positive and negative $\pi$

## Avoiding the square root

I'll look back at Rodrigues' rotation formula now, and write it in terms of the *half-angle* rather than the full angle $\phi = |\omega| \Delta{t}$. This will turn out helpful later, trust me.

$$
\begin{align}
R(\omega\Delta{t}) &= \mathbb{I} + \sin(\phi) [\hat{n}]_\times + (1 - \cos(\phi)) [\hat{n}]_\times^2 \\
&= \mathbb{I} + 2 \cos(\phi/2) \sin(\phi/2) [\hat{n}]_\times + (1 - (\cos^2(\phi/2)-\sin^2(\phi/2))) [\hat{n}]_\times^2 \\
&= \mathbb{I} + 2 \cos(\phi/2) \sin(\phi/2) [\hat{n}]_\times + 2\sin^2(\phi/2) [\hat{n}]_\times^2 \\
&= \mathbb{I} + 2 \cos^2(\phi/2) \left[\tan(\phi/2) [\hat{n}]_\times + \tan^2(\phi/2) [\hat{n}]_\times^2\right] \\
&= \mathbb{I} + 2 \frac{\tan(\phi/2) [\hat{n}]_\times + \tan^2(\phi/2) [\hat{n}]_\times^2}{1 + \tan^2(\phi/2)} \\
\end{align}
$$
Now I will intentionally apply the same approximation as for the quaternion case above, that is $\phi = 2 \arctan(|\omega|\Delta{t}/2)$, and so $\tan(\phi/2) = |\omega|\Delta{t}/2$:
$$
R(\omega\Delta{t}) \approx \mathbb{I} + 2 \frac{[\frac{\omega\Delta{t}}{2}]_\times + [\frac{\omega\Delta{t}}{2}]^2_\times}{1 + |\frac{\omega\Delta{t}}{2}|^2}
$$
and so:
$$
R(\omega\Delta{t}) \cdot v \approx v + 2 \frac{\frac{\Delta{t}}{2} (\omega \times v) + \frac{\Delta{t}}{2}^2 \omega \times (\omega \times v)}{1 + |\frac{\omega\Delta{t}}{2}|^2}
$$
And there we have it: an approximate rotation of a vector by a given angular velocity and timestep, without needing any non-elementary arithmetic operations! While the angle turned through is only correct to second-order:
$$
\Theta = |\omega| \Delta{t} + \frac{|\omega|^3\Delta{t}^3}{12} + O(|\omega|^5\Delta{t}^5))
$$
the axis is the same as the exact rotation, and the length-preserving property of rotation operators is reproduced without requiring a square root evaluation to compute the length of an intermediate result. 

Now, this is all good for operating on individual vectors, but what if we want to update an entire orientation? That ideally needs a quaternion to multiply by, not a rotation matrix. I'll try applying the same half-angle trick to the expression for a quaternion:
$$
\begin{align}
q(\omega \Delta{t}) &= \left[\cos(|\omega|\Delta{t}/2), \hat{n}\sin(|\omega|\Delta{t}/2)\right] \\
&= \left[\cos^2(|\omega|\Delta{t}/4)-\sin^2(|\omega|\Delta{t}/4), 2 \hat{n}\sin(|\omega|\Delta{t}/4)\cos(|\omega|\Delta{t}/4)\right] \\
&= \cos^2(|\omega|\Delta{t}/4) \left[1 -\tan^2(|\omega|\Delta{t}/4), 2 \hat{n}\tan(|\omega|\Delta{t}/4)\right] \\
&= \frac{\left[1 -\tan^2(|\omega|\Delta{t}/4), 2 \hat{n}\tan(|\omega|\Delta{t}/4)\right]}{1 + \tan^2(|\omega|\Delta{t}/4)} \\
&\approx \frac{\left[1 - \left(\frac{|\omega|\Delta{t}}{4}\right)^2, 2 \frac{\omega\Delta{t}}{4}\right]}{1 + \left(\frac{\omega|\Delta{t}}{4}\right)^2}
\end{align}
$$
This time the approximation is $\Theta = 4 \arctan(|\omega|\Delta{t}/4)$, which is also second-order accurate but with a third-order coefficient 4 times smaller:
$$
\Theta = |\omega| \Delta{t} + \frac{|\omega|^3\Delta{t}^3}{48} + O(|\omega|^5\Delta{t}^5))
$$
This would be the form I would recommend for updating quaternion orientations if speed is of the essence, with the above form used for rotating vectors. 

## Benchmarking

In some basic testing I have carried out, the methods I have suggested are as fast or faster than others, and can be substantially faster than constructing a full quaternion.

For rotating a vector, doing so by computing the full quaternion and performing a quaternion-vector conjugation took 28.9 nanoseconds. Computing my recommended approximate quaternion and conjugating took 22.3 ns. Doing the "first-order and normalise" took 20.1 ns and my recommended approximate vector rotation took 18.4 ns. My suggested method is both the fastest, and has accuracy advantages (i.e. it reproduces the correct axis) over the next fastest method, the "first-order and normalise". It's over a third faster than the full calculation, so is a good option where rotation rates are relatively low or execution speed is the priority over maximum accuracy.

For computing an orientation quaternion (not including multiplying with another quaternion, as this would be the same for all of the methods), the "exact" method took 21.8 ns, "first-order and normalise" took 16.7 ns and my approximation took 15.0 ns. Here, the "first-order and normalise" method *does* reproduce the correct axis, but its angle error is slightly worse than that of my method, as well as being marginally slower.