When simulating rigid-body dynamics in three dimensions, we quickly run into an issue that's not present in two: rotating objects in three dimensions behave in ways that can be surprising and are kind of weird. In this post I'll present a theoretical framework for describing this behaviour, then look at some existing approaches to simulating the motion of rotating rigid bodies and finally suggest a method with improved energy conservation properties, which I believe to be novel.

## First Principles

In this post I'll describe the three classes of rotating objects, how they behave and how we can treat them. The dynamics of a freely rotating body can be derived from the following:
1. An object has a *moment of inertia* $I$, the rotational analog of mass, which is constant and fixed to the object itself. The moment of inertia is a $3 \times 3$ symmetric positive-definite matrix, meaning that in some frame of reference it is a diagonal matrix with positive entries. This frame of reference rotates with the body, so I'll denote it the *body frame*. The coordinate axes of the body frame (the eigenvalues of the moment of inertia) are called the body's principal axes. It's not really relevant here, but it makes most sense to place the origin of the body frame at the object's centre of mass, as this is the point a free object will rotate around. As it is symmetric positive-definite, The moment of inertia matrix is always invertible. I will use its inverse in this post a lot, so I will give it its own symbol: $J = I^{-1}$.
2. An object has an *angular momentum* $L$, which is a conserved quantity for free motion (torque-free, drag-free, contact-free, you know the drill). Its length is constant in any frame of reference, and its direction is constant in any frame of reference that is not itself rotating (being a physicist, I may call such a frame of reference the *lab frame*, and I'll use this fairly interchangeably with *world-space*).
3. The object's orientation $R$, at any instantaneous point in time, rotates with an angular velocity of $\omega = J L$. In the lab frame, $L$ is constant but $I$ (and therefore $J$) is not. In the body frame, $I$ and $J$ are constant but $L$ is not. $\omega$ is not necessarily constant in either frame of reference, though it can be if $L$ is aligned with one of the object's principal axes. 
4. The final conserved quantity is the rotational energy, half the dot product of the angular velocity and the angular momentum: $E = \frac{1}{2} \omega \cdot L = \frac{1}{2} LJL = \frac{1}{2} \omega I \omega$.

### Rotations 

Something I'll use throughout this post is Rodrigues' rotation formula, which expresses a rotation matrix as a function of its "scaled axis" representation (i.e. $\vec{\phi} = \omega dt = |\phi| \hat{n}$). In several equivalent forms:
$$
\begin{align*}
R(\vec{\phi}) &= \mathbb{1} + (\sin{|\phi|}) [\hat{n}]_\times + (1 - \cos{|\phi|}) [\hat{n}]_\times^2 \\
&= \mathbb{1} + (\mathrm{sinc}{|\phi|}) [\vec{\phi}]_\times + \frac{1}{2}\mathrm{sinc}^2{\frac{|\phi|}{2}} [\vec{\phi}]_\times^2 \\
&= \mathbb{1} + \left(\mathrm{sinc}{\frac{|\phi|}{2}}\cos{\frac{|\phi|}{2}}\right) [\vec{\phi}]_\times + \frac{1}{2}\mathrm{sinc}^2{\frac{|\phi|}{2}} [\vec{\phi}]_\times^2
\end{align*}
$$
A first interesting result here is that, if you're happy to conflate $\phi/2$ with $\tan(\phi/2)$ (which is common for many finite-difference rotation methods that do first order update followed by normalisation), you can get to a length-preserving rotation operation without the need for any transcendental (e.g. trigonometric or exponential) or even square root calculations! (for more on this see [[Fast, Approximate Rotation Operators]])
$$
\begin{align*}
\cos{x} &= \frac{1}{\sqrt{1 + \tan^2{x}}} \\
\sin{x} &= \frac{\tan{x}}{\sqrt{1 + \tan^2{x}}} \\
R(\vec{\phi}) &= \mathbb{1} + 2\frac{\tan\frac{|\phi|}{2}}{1 + \tan^2\frac{|\phi|}{2}} \frac{[\vec{\phi}]_\times}{|\phi|} + 2 \frac{\tan^2\frac{|\phi|}{2}}{1 + \tan^2\frac{|\phi|}{2}}\frac{[\vec{\phi}]_\times^2}{|\phi|^2} \\
&\approx \mathbb{1} + 2\frac{\frac{|\phi|}{2}}{1 + \left(\frac{|\phi|}{2}\right)^2} \frac{[\vec{\phi}]_\times}{|\phi|} + 2\frac{\frac{|\phi|}{2}^2}{1 + \left(\frac{|\phi|}{2}\right)^2}\frac{[\vec{\phi}]_\times^2}{|\phi|^2} \\
&= \mathbb{1} + \frac{[\vec{\phi}]_\times + \frac{1}{2} [\vec{\phi}]_\times^2}{1 + \frac{\phi^2}{4}}
\end{align*}
$$
For a typical "rotate vector according to a given angular velocity" situation we have:
$$
R(\vec\omega \Delta t) \cdot \vec{v} = \vec{v} + 2\frac{\vec{\omega}\times\vec{v}\ \frac{\Delta{t}}{2} + \vec{\omega}\times\left(\vec{\omega} \times \vec{v}\right)\left(\frac{\Delta{t}}{2}\right)^2 }{1 + \omega^2 (\Delta{t} / 2)^2}
$$
That's two cross products, a squared magnitude, some additions and a vector-scalar division, for a length-preserving rotation operation! The only compromise is that it slightly distorts time, rotating for a modified time $\omega \Delta t' = 2 \arctan(\omega \Delta{t}/2)$, meaning about 5% of the rotation angle is lost if you rotate by an eighth of a full turn in a single iteration, or about 15% for a quarter-turn per iteration.

As noted earlier this is an equivalent loss of angular velocity to the explicit integration of rotations by doing $q(t + \Delta{t}) = q(t) + 1/2 \dot{q} \Delta{t}$ and then renormalising the quaternion.

Something that's rotating _very_ fast will be limited to a maximum of a half-turn per iteration, but that's already the case for many low-order rotation approximations. I assume this formula isn't novel, but I've never seen it before and I was surprised to see it drop out like that! I'd always assumed either trig functions or square roots would be required for length-preserving rotations.

_Anyway_, back to rigid-body dynamics. I'll start by going over some background, deriving the usual formulas from the conservation laws I proposed above.

## Euler's Equations

Call the (constant) body-frame inertia and its inverse $I_b$ and $J_b$ respectively. If the body has an orientation (or "attitude") given by the rotation matrix $R$, then in the lab frame, we can compute the world-space inertia matrix (and its inverse):
$$
\begin{align*}
I_w = R\ I_b\ R^T \\
J_w = R\ J_b\ R^T
\end{align*}
$$
$R$ is (by definition) not a constant for a rotating body, so the lab-frame inertias may not be either.
We can also calculate the body-space angular momentum:
$$
L_b = R^T L_w
$$
The evolution of the angular velocity follows from angular momentum conservation :
$$
\begin{align*}
\frac{d\omega_w}{dt} &= \frac{d(J_wL_w)}{dt} \\
&= \frac{dJ_w}{dt} L_w + J_w \frac{dL_w}{dt} \\
&= \frac{dJ_w}{dt} L_w
\end{align*}
$$
$$
\begin{align*}
\frac{dJ_w}{dt} &= \frac{dR}{dt} J_b R^T + R J_b \frac{dR}{dt}^T \\
&= [\omega_w]_\times R J_b R^T - R J_b R^T [\omega_w]_\times
\end{align*}
$$
$$
\begin{align*}
\frac{d\omega_w}{dt} = \frac{dJ_w}{dt} L_w &= \omega_w \times (R J_b R^T L_w) - R J_b R^T (\omega_w \times L_w) \\
&= \omega_w \times \omega_w - J_w (\omega_w \times L_w) \\
&= - J_w (\omega_w \times L_w)
\end{align*}
$$
In body coordinates, similarly:
$$
\begin{align*}
\frac{d\omega_b}{dt} = J_b \frac{dL_b}{dt} = J_b \frac{dR^T}{dt}L_w = -J_b R^T (\omega_w \times L_w) = -J_b (\omega_b \times L_b)
\end{align*}
$$
In the above I have used the derivative of the body's orientation, referencing point 3 above and Rodrigues' formula:
$$
\frac{dR}{dt} = \lim_{\Delta t \rightarrow 0} \left[
\frac{R(t + \Delta t) - R(t)}{\Delta t} 
= \frac{Q(\omega \Delta t) R(t) - R(t)}{\Delta t}
= \frac{Q(\omega \Delta t) - \mathbb{1}}{\Delta t} R(t)
\right] =  [\omega]_\times R(t)
$$
$$
\frac{dR^T}{dt} = \lim_{\Delta t \rightarrow 0} \left[
\frac{R^T(t + \Delta t) - R^T(t)}{\Delta t} 
= \frac{R(t)^T Q^T(\omega \Delta t) - R^T(t)}{\Delta t}
= R(t)^T\frac{ Q(-\omega \Delta t) - \mathbb{1}}{\Delta t}
\right] =  -R^T(t) [\omega]_\times 
$$

The differential equation for the angular velocity comes out the same in both cases:
$$ \frac{d\omega}{dt} = - J (\omega \times L) $$
These are the famous Euler's equations for rigid body dynamics (not to be confused with _Euler's equation_, which describes incompressible fluids. Give someone else a go, Leonhard). In their usual formulation they're stated component-wise, in the body frame:
$$
\begin{align*}
I_1 \dot{\omega}_1 + (I_3 - I_2) \omega_2 \omega_3 = 0\\ 
I_2 \dot{\omega}_2 + (I_1 - I_3) \omega_3 \omega_1 = 0\\ 
I_3 \dot{\omega}_3 + (I_2 - I_1) \omega_1 \omega_2 = 0
\end{align*}
$$
In the body frame, the constancy of the moment of inertia lets us get a differential equation for the evolution of the angular momentum (in world space, this is just the conservation law $dL = 0$):
$$
\frac{dL_b}{dt} = \frac{dI_b\omega_b}{dt} = I_b \frac{d\omega_b}{dt} = - I_b J_b (\omega_b \times L_b) = - \omega_b \times L_b = - (J_b L_b) \times L_b
$$
We can write an analogous set of equations for the angular momentum in body space:
$$
\begin{align*}
\dot{L}_1 + (J_2 - J_3) L_2 L_3 = 0\\ 
\dot{L}_2 + (J_3 - J_1) L_3 L_1 = 0\\ 
\dot{L}_3 + (J_1 - J_2) L_1 L_2 = 0
\end{align*}
$$
I like these equations because they involve only quantities that are constant in one frame or the other: as $L_b$ changes in the body frame, it's clear that the body must have rotated in world space such that $L_w$ remains unchanged. 

Note that just knowing the evolution of the angular momentum (or angular velocity) vector in body space does not _fully_ specify the solution we're looking for, that is the evolution of the object's orientation in world space. The body-space solution is blind to rotations of the body around the angular momentum, as these do not affect the body-space angular momentum.

A full solution can be obtained by choosing a unit vector orthogonal to the angular momentum, call it $v$. Then, the motion of the body consists of choosing the rotation which both fulfils the Euler equations *and* maintains the direction of $v$ in world-space. The evolution of $v$ in world-space is constrained to the plane perpendicular to the angular momentum, so it can be described as a rotation around the angular momentum. Composing the rotation induced by evolving the Euler equations (to keep $L$ and $v$ constant in world-space) with a rotation around $L$ itself (to match the needed evolution of $v$ in world-space) produces the full dynamics.

A full solution is treated in the excellent paper [Numerical implementation of the exact dynamics of
free rigid bodies](https://arxiv.org/pdf/cond-mat/0607529). The authors advocate for implementing the exact solutions in terms of Jacobi functions, which I will not advocate for here. I might do some noodling with them at some point later though. 

Now I've covered the general shape of the problem, I will focus on the three particular special cases: the spherical top, the symmetric top and the asymmetric top.

## The Spherical Top

If we're just thinking about the most symmetric objects, everything is fine. Spheres and cubes fit into this category, as do other platonic solids. The moment of inertia of an object in this category is isotropic - a scalar multiple of the identity matrix - which means that their angular velocity and angular momentum are always parallel to one another. In the absence of external torques, these objects will just spin around that common axis forever perfectly happily, with constant angular velocity.

Euler's equations become trivial in this case, because $J_1 = J_2 = J_3$:
$$
\dot{L}_1 = 0 \qquad \dot{L_2} = 0 \qquad \dot{L_3} = 0
$$
To advance a simulated body by a time $\Delta t$, we can just rotate it by $\omega \Delta{t}$:
$$
R(t + \Delta t) = Q(\omega \Delta t) R(t) = Q(J L \Delta t) R(t)
$$

## The Symmetric Top

The next simplest category of objects are those with two of their principal angular inertia values identical, and one different, exemplified by the cylinder. You could imagine that it makes a difference whether the cylinder is a long rod (with its unique principal inertia greater than the others) or a disc (with unique inertia less than the other two), and while these cases look quite different the mathematics works out mostly the same in both cases.
The cylindrical object will, if its angular momentum is oblique (neither parallel nor perpendicular) to its special axis, spin in a precessing pattern. It's easiest to think about this motion as spinning at a particular speed around the axis of the object, composed with a rotation of that axis around the angular momentum vector (which, as a conserved quantity, remains constant).

Looking at the Euler equation, with $J_1 = J_2 = J_\perp$ and $J_3 = J_\parallel$ :
$$
\begin{align*}
\dot{L}_1 + (J_2 - J_3) L_2 L_3 &= 0 \\
\dot{L}_2 + (J_3 - J_1) L_3 L_1 &= 0 \\
\dot{L}_3 + (J_1 - J_2) L_1 L_2 &= 0 
\end{align*}
$$
$$
\begin{align*}
\dot{L}_1 - (J_\parallel - J_\perp) L_3 L_2 &= 0 \\
\dot{L}_2 + (J_\parallel - J_\perp) L_3 L_1 &= 0 \\
\dot{L}_3 &= 0 
\end{align*}
$$
The angular momentum component along the special axis, $L_3$, is constant. 
The other two components look like a planar oscillation with frequency $\Omega = (J_\parallel - J_\perp) L_3$:
$$
\begin{align*}
\dot{L}_1 &= \Omega L_2 \\
\dot{L}_2 &= - \Omega L_1
\end{align*}
$$
Substituting each equation into the time derivative of the other yields an oscillator ODE for each component, indicating that the "perpendicular" component of the angular momentum, in body space, rotates around its special axis with frequency $\Omega$.
$$
\begin{align*}
\ddot{L}_1 &= \Omega \dot{L}_2  = -\Omega^2 L_1\\
\ddot{L}_2 &= - \Omega \dot{L}_1 = -\Omega^2 L_2
\end{align*}
$$
To advance the simulated body in world space, rotate by $\Omega \Delta t$ about the object's unique axis ($R(t) \hat{z}$ in world space). I note here that this axial angular velocity is that obtained from the world-space angular momentum and a modified inverse inertia tensor $(J_w - J_\perp)$ - where the subtraction of a scalar from a matrix implicitly multiplies the scalar by the identity matrix, such that $(J_w - J_\perp) L_w = J_w L_w - J_\perp L_w$. 
$$\Omega R(t) \hat{z} = (J_\parallel - J_\perp) L_3  R(t) \hat{z} = (J_w - J_\perp) L_w$$
Then rotate the body around its angular momentum vector by an angle $J_\perp L_w \Delta t$. Rotating the body around the angular momentum in world space leaves the body-space angular momentum unchanged, so we still satisfy Euler's equations. Note that the "precession" angular velocity $J_\perp L_w$ and the previous, "axial" angular velocity $(J_w - J_\perp)L_w$ sum to $\omega_w = J_w L_w$. That means the instantaneous rotation of the body has angular velocity $\omega_w$, as required by point 3 above.

All in all, we have:
$$
R(t + \Delta t) = Q(J_\perp L_w \Delta t) Q((J_w - J_\perp)L_w \Delta t) R(t)
$$
This amounts to a "splitting" the contribution of the total angular velocity $J_w L_w$ into the "axial" rotation $Q[(J-J_\perp) L \Delta t]$ (which solves the Euler equation, rotating the angular momentum in body space (or equivalently the body in world space) in a way that preserves the rotational energy of the body) and the "precessional" rotation $Q[J_\perp L \Delta t]$, which rotates the body around the angular momentum vector and therefore does not affect any of our conserved quantities. 

## The Asymmetric Top

The most complicated, and most general type of body is referred to as the asymmetric top. Bodies in this category have three distinct principal moments of inertia, and can "tumble" in ways that appear unpredictable and chaotic. If you're reading this, and have got this far, you're probably already aware of the [Dzhanibekov effect/Tennis racquet theorem/intermediate axis theorem](https://en.wikipedia.org/wiki/Tennis_racket_theorem). If not, the linked Wikipedia page has a bunch of nice demonstrations of it. 

The long and short of it is: an asymmetric object that's rotating with angular momentum mostly aligned to its "shortest" or "longest" axis will precess in a little ellipse, similar to the way a symmetric object precesses in a circle. As the angular momentum approaches the intermediate axis, the nice ellipse distorts into a weird shape (Wikipedia compares it to a taco). If its angular momentum is close enough to the "intermediate" axis, the precessions are less like a steady, gentle rocking back and forth, and more like a periodic, violent inversion of the object. An analytic solution for this motion (the equivalent to what I presented above for the more symmetric classes of body) _does exist_, but it uses significantly more exotic special functions than the trig functions I've been relying on thus far (in Rodrigues' formulas for rotation operators). The paper [Numerical implementation of the exact dynamics of free rigid bodies](https://arxiv.org/abs/cond-mat/0607529) goes through the exact solution using the Jacobi elliptic functions $\mathrm{cn}$, $\mathrm{sn}$ and $\mathrm{dn}$ which roughly resemble $\cos$, $\sin$ and $1$ respectively. Their behaviour is modified by a parameter $0 \leq m \leq 1$ which quantifies how "taco-like" the precessions are, and increases as the angular momentum gets closer to the intermediate axis. 

I'll review some prior art for solving this general case, talking through how well they obey the conservation laws I wrote at the top of this post.

### Forward Euler

The most naïve approach possible is good old Forward Euler, no mitigations. 
$$
\begin{align*}
\omega_1 &= \omega_0 - J(\omega_0 \times L_0) \Delta{t} \\
R_1 &= Q(\omega_1 \Delta{t}) R_0
\end{align*}
$$
This was approximately the method used in the earliest versions of Avian physics (back when it was `bevy_xpbd`, and up to Avian 0.1.0) It was abandoned on account of being massively unstable for asymmetric objects: look at the evolution of the magnitude of the angular momentum (a quantity that should be conserved!):
$$
\begin{align*}
L_1^2 = |I\omega_1|^2 &= \omega_0 I^2 \omega_0 - 2 \omega_0 I^2 J (\omega_0 \times L_0) \Delta{t} + (J(\omega_0 \times L_0))I^2(J(\omega_0 \times L_0)) \Delta{t}^2 \\ 
&= L_0 \cdot L_0 - 2 L_0 \cdot (\omega_0 \times L_0) \Delta{t} + |\omega_0 \times L_0|^2 \Delta{t}^2 \\
&= L_0^2 + k^2 \Delta{t}^2
\end{align*}
$$
The quantity $\omega_0 \times L_0$ is going to appear repeatedly, so I've given it a shorthand name, $k$. As the squared magnitude of a vector, $k^2$ is non-negative, so the magnitude of the angular momentum increases monotonically. This is actually worse than it looks: As $k = (JL) \times L$ is quadratic in the magnitude of the angular momentum, $L_0^2 \approx L_0^2 + L_0^4 \Delta{t}^2$ - the growth in the angular momentum is _faster than exponential_.

Needless to say, energy conservation also doesn't come out of this well:
$$
\begin{align*}
2E_1 = \omega_1 \cdot L_1 &= \omega_0 I \omega_0 - 2 \omega_0 I J (\omega_0 \times L_0) \Delta{t} + (J(\omega_0 \times L_0))I(J(\omega_0 \times L_0)) \Delta{t}^2 \\
&= \omega_0 \cdot L_0 - 2 \omega_0 \cdot (\omega_0 \times L_0) \Delta{t} + (\omega_0 \times L_0)J(\omega_0 \times L_0) \Delta{t}^2 \\
&= 2E_0 + kJk \Delta{t}^2
\end{align*}
$$
Its not as obvious, but because $J$ is positive-definite, $kJk$ is positive whenever $k$ is nonzero. The energy also blows up to infinity. Not ideal.

### The "Jolt" method

After some time using an implicit method for its rotational integration (which I won't go into here - being implicit makes it harder to analyse!) in response to the instability of the explicit method, Avian has now switched to using a method borrowed from the Jolt physics engine. This method is conceptually pretty simple: just apply an explicit update to the angular momentum, then renormalise it!
$$
\begin{align*}
L' &= L_0 - (\omega_0 \times L_0)\ \Delta{t} \\
L_1 &= L' \sqrt{\frac{L_0^2}{L'^2}} = \frac{L'}{\sqrt{1 + \Delta{t}^2 k^2/L_0^2}} \\
\omega_1 &= J L_1 \\
R_1 &= Q(\omega_1 \Delta{t}) R_0 
\end{align*}
$$
The magnitude of the angular momentum is explicitly conserved, so that resolves the first of the blow-ups we saw in the explicit method. Looking at the energy:
$$
\begin{align*}
2E_1 = \omega_1 \cdot L_1 &= \frac{\left[L_0 - k \Delta{t}\right] J \left[L_0 - k \Delta{t}\right]}{1 + \Delta{t}^2 k^2/L_0^2} \\
&= \frac{L_0 J L_0 + kJk \Delta{t}^2}{1 + \Delta{t}^2 k^2/L_0^2} \\
&= \frac{2E_0 + kJk \Delta{t}^2}{1 + \Delta{t}^2 k^2/L_0^2} \\
&= 2E_0 + \frac{kJk - 2 E_0 k^2/L_0^2}{1 + \Delta{t}^2 k^2/L_0^2} \Delta{t}^2 \\
2 \Delta{E} &= (J_k - J_L) L_0^2 \frac{\Delta{t}^2 k^2/L_0^2}{1 + \Delta{t}^2 k^2/L_0^2}
\end{align*}
$$
Here I have introduced the notation $kJk = J_k k^2$ and similarly $L_0 J L_0 = J_L L_0^2$. The $J_k$ and $J_L$ values are the effective scalar magnitude of the matrix $J$, in the direction of $k$ and $L_0$ respectively, which are necessarily bounded by the eigenvalues of $J$ and will generally be some weighted average of those eigenvalues. The timestep-dependent factor $\frac{\Delta{t}^2 k^2/L_0^2}{1 + \Delta{t}^2 k^2/L_0^2}$ looks like a quadratic for small values of its parameter, $|k|\Delta{t}/|L_0| = |\omega_0| \Delta{t} \sin{\theta}$ but saturates to 1 as the parameter becomes large. Importantly, the sign of $\Delta{E}$ is now not necessarily positive. $(J_k - J_L)$ can be either positive or negative, and qualitative experimentation indicates that the overall effect over several timesteps is a reasonably good conservation of energy. In particular, because $J_\mathrm{min} \leq J_L \leq J_\mathrm{max}$ and $2E = J_L L^2$ and because this method explicitly conserves $L^2$, the energy is bounded from above and below: $\frac{1}{2} J_\mathrm{min} L_0^2 \leq E \leq \frac{1}{2} J_\mathrm{max} L_0^2$. I speculate that there's a stabilising effect from the $(J_k - J_L)$ term, but I've not proven this rigorously, or even tested it particularly thoroughly. 

### Rotating the Angular Momentum

Looking at the Jolt method, the initial modification of the angular momentum looks a bit like it's a first-order approximation to rotating $L_0$ by $-\omega_0 \Delta{t}$. So, how about doing this:
$$
\begin{align*}
L_1 &= Q(-\omega_0 \Delta{t}) L_0 \\
\omega_1 &= J L_1 \\
R_1 &= Q(\omega_1 \Delta{t}) R_0
\end{align*}
$$
Once again our angular momentum's magnitude is conserved, so let's look straight at the energy. I'll employ Rodrigues' rotation formula, with the shorthand $s = \mathrm{sinc}(\omega_0 t / 2)$ and $c = \cos(\omega_0 t / 2)$:
$$
\begin{align*}
L_1 &= L_0 - s c (\omega_0 \times L_0) \Delta{t} + \frac{1}{2} s^2 \omega_0 \times (\omega_0 \times L_0) \Delta{t}^2 \\
&= L_0 - s c k \Delta{t} + \frac{1}{2} s^2 (\omega_0 \times k) \Delta{t}^2
\end{align*}
$$
Using the following expansions:
$$
\begin{align*}
(\omega_0 \times k) &= (\omega_0 \cdot L_0) \omega_0 - \omega_0^2 L_0 \\
k \cdot J(\omega_0 \times k) &= (\omega_0 \cdot L_0) (\omega_0Jk) \\
(\omega_0 \times k)\cdot J(\omega_0 \times k) &= [(\omega_0 \cdot L_0) \omega_0 - \omega_0^2 L_0] J [(\omega_0 \cdot L_0) \omega_0 - \omega_0^2 L_0] \\ 
&= (\omega_0\cdot L_0)^2 (\omega_0 J \omega_0) - (\omega_0\cdot L_0) \omega_0^4
\end{align*}
$$
I get to:
$$
\begin{align*}
2 E_1 = L_1 J L_1 = L_0 J L_0 &- 2 s c (JL_0)\cdot k\Delta{t} \\
&+ s^2 (JL_0) \cdot (\omega_0 \times k) \Delta{t}^2 \\
&+ s^2 c^2 (kJk) \Delta{t}^2 \\
&- s^3 c (Jk)\cdot(\omega_0\times k) \\
&+ \frac{1}{4} s^4 (\omega_0 \times k)J(\omega_0 \times k) \Delta{t}^4 \\
= 2 E_0 &+ s^2 c^2 (kJk) \Delta{t}^2 \\
&- s^3 c (\omega_0 \cdot L_0) (\omega_0 J k) \Delta{t}^3 \\
& + \frac{1}{4} s^4 \left[(\omega_0\cdot L_0)^2 (\omega_0 J \omega_0) - (\omega_0\cdot L_0) \omega_0^4\right] \Delta{t}^4
\end{align*}
$$

To lowest order in time,
$$\Delta{E} = \frac{1}{2}\mathrm{sinc}^2\left(\frac{\omega_0\Delta{t}}{2}\right)\cos^2\left(\frac{\omega_0\Delta{t}}{2}\right) kJk \Delta{t}^2 = \frac{1}{2} kJk \Delta{t}^2.$$
We still have the monotonic growth in energy shown by the explicit method, not the $(J_k - J_L)$ factor we had for the Jolt method (although now that $L^2$ is conserved, at least the maximum energy is bounded).

This isn't what we wanted, and implies the Jolt method is actually doing something different to just "low order approximation to rotating around $\omega_0$". It took me a long time and a lot of scribbling to figure out what the difference is, and the key is that the normalisation step, $L_1 = L' \sqrt{L_0^2 / L'^2}$, shrinks *all* components of $L'$ isotropically, whereas a *true* rotation leaves the component of a vector parallel to the axis unchanged. Then I realised that while $L' = L_0 - \omega_0 \times L_0 \Delta{t}$ is *indeed* a first-order approximation to a rotation around $\omega_0$, it's not a one-to-one correspondence: treated as a linear operator, the cross product with $L_0$ has a zero eigenvalue (with the corresponding eigenvector being $L_0$ itself). That means that for *any* real scalar value $J'$, the above expression is a first-order approximation to a rotation around $\omega_0 - J' L_0$, with the same $\Delta{t}$:
$$
\begin{align*}
L' &= L_0 - (\omega_0 - J' L_0) \times L_0 \Delta{t}\\
&= L_0 - \omega_0 \times L_0 \Delta{t} + J' (L_0 \times L_0) \Delta{t}\\
&= L_0 - \omega_0 \times L_0 \Delta{t}
\end{align*}
$$
### The Jolt Method in Rotation Form (Properly)

All that's now necessary to reproduce the Jolt method is to figure out the right value of $J'$. Then the Jolt method is approximately a rotation with angular velocity $\Omega = \omega_0 - J' L_0$. The key insight here is that if $L_1 = L' \sqrt{L_0^2 / L'^2}$, shrinks all of the components of $L'$, but the component parallel to the axis must be unchanged, then that component must be zero!
$$
\begin{align*}
L' \cdot \Omega &= 0 \\
(L_0 - \Omega \times L_0) \cdot \Omega &= 0 \\
L_0 \cdot \Omega &= 0 \\
L_0 \cdot (\omega_0 - J' L_0) &= 0 \\
J' &= \frac{\omega_0 \cdot L_0}{L_0^2} = \frac{L_0 J L_0}{L_0^2} = J_L
\end{align*}
$$
So we have $\Omega = \omega_0 - J_L L_0 = (J - J_L) L_0$. Cool! What's our energy looking like now? First, angular momentum:
$$
\begin{align*}
L_1 &= L_0 - s c (\Omega \times L_0) \Delta{t} + \frac{1}{2} s^2 \Omega \times (\Omega \times L_0) \Delta{t}^2 \\
&= L_0 - s c (\omega_0 \times L_0) \Delta{t} + \frac{1}{2} s^2 \Omega \times (\omega_0 \times L_0) \Delta{t}^2 \\
&= L_0 - s c k \Delta{t} + \frac{1}{2} s^2 (\Omega \times k) \Delta{t}^2
\end{align*}
$$
Where I've used the fact that $\Omega \times L_0 = \omega_0 \times L_0 = k$ to simplify. Now the energy:
$$
\begin{align*}
2E_1 &= L_1JL_1 \\
&= L_0 J L_0 + s^2 c^2 kJk \Delta{t}^2 + s^2 L_0 J (\Omega \times k) \Delta{t}^2 - s^3 c k J (\Omega\times k) \Delta{t}^3 + \frac{1}{4} s^4 (\Omega \times k)J(\Omega \times k) \Delta{t}^4 \\
&= 2E_0 + s^2 (c^2 kJk - J' \omega_0 \cdot(L_0 \times k)) \Delta{t}^2 - s^3 c (\Omega \cdot L) (\omega_0 J k) \Delta{t}^3 + \frac{1}{4}s^4 (\Omega \times k)J(\Omega \times k) \Delta{t}^4 \\
&= 2E_0 + s^2 (c^2 J_k - J') k^2 \Delta{t}^2 - s^3 c (\Omega \cdot L) (\omega_0 J k) \Delta{t}^3 + \frac{1}{4}s^4 (\Omega \times k)J(\Omega \times k) \Delta{t}^4 \\
&= 2E_0 + s^2 (c^2 J_k - J') k^2 \Delta{t}^2 - s^3 c (\Omega \cdot L) (\omega_0 J k) \Delta{t}^3 + \frac{1}{4}s^4 \left[(\Omega \cdot L) \omega_0 - (\Omega \cdot \omega_0) L\right]J\left[(\Omega \cdot L) \omega_0 - (\Omega \cdot \omega_0) L\right] \Delta{t}^4 \\
\end{align*}
$$
For the $\Delta{t}^3$ term, I have used the following simplification:
$$
k \times Jk = (Jk) \times (L \times \omega) = (\omega Jk) L - (\omega \cdot k) \omega = (\omega J k) L
$$
Substituting the value $J' = J_L$, and therefore that $(\Omega \cdot L) = 0$ we have:
$$
\begin{align*}
2 \Delta{E} &= s^2 (c^2 J_k - J_L) k^2 \Delta{t}^2 + \frac{1}{4}s^4 (\Omega \cdot \omega_0)^2 LJL \Delta{t}^4 \\
&= s^2 (c^2 J_k - J_L) k^2 \Delta{t}^2 + \frac{1}{4}s^4 (\Omega \cdot \omega_0)^2 (\omega_0 \cdot L_0) \Delta{t}^4 \\
&= s^2 (c^2 J_k - J_L) k^2 \Delta{t}^2 + \frac{1}{4}s^4 \left[\Omega^2 \omega_0^2 - |\Omega\times\omega_0|^2\right] (\omega_0 \cdot L_0) \Delta{t}^4 \\
&= s^2 (c^2 J_k - J_L) k^2 \Delta{t}^2 + \frac{1}{4}s^4 \left[\Omega^2 \omega_0^2 - J_L^2 k^2\right] (\omega_0 \cdot L_0) \Delta{t}^4 \\
&= s^2 (c^2 J_k - J_L) k^2 \Delta{t}^2 + \frac{1}{4}s^4\Omega^2 \omega_0^2 J_L L_0^2 \Delta{t}^4 - \frac{1}{4}s^4  J_L^2 k^2 (\omega_0 \cdot L_0) \Delta{t}^4  \\
&= s^2 (c^2 J_k - J_L) k^2 \Delta{t}^2 + \frac{\Omega^2\Delta{t}^2}{4}s^2 s^2 \omega_0^2 J_L L_0^2 \Delta{t}^2 - \frac{1}{4}s^4  J_L^2 k^2 (\omega_0 \cdot L_0) \Delta{t}^4  \\
&= s^2 (c^2 J_k - J_L) k^2 \Delta{t}^2 + \sin^2\left(\frac{\Omega \Delta{t}}{2}\right) s^2  J_L L_0^2\omega_0^2 \Delta{t}^2 - \frac{1}{4}s^4  J_L^2 k^2 (\omega_0 \cdot L_0) \Delta{t}^4 \\
&= s^2 (c^2 J_k - J_L) k^2 \Delta{t}^2 + \sin^2\left(\frac{\Omega \Delta{t}}{2}\right) s^2  J_L [k^2 + (\omega_0 \cdot L_0)^2] \Delta{t}^2 - \frac{1}{4}s^4  J_L^2 k^2 (\omega_0 \cdot L_0) \Delta{t}^4 \\
&= s^2 c^2 (J_k - J_L) k^2 \Delta{t}^2 + \sin^2\left(\frac{\Omega \Delta{t}}{2}\right) s^2  J_L (\omega_0 \cdot L_0)^2 \Delta{t}^2 - \frac{1}{4}s^4  J_L^2 k^2 (\omega_0 \cdot L_0) \Delta{t}^4 \\
&= s^2 c^2 (J_k - J_L) k^2 \Delta{t}^2 + \frac{1}{4} s^4  \Omega^2 J_L (\omega_0 \cdot L_0)^2 \Delta{t}^4 - \frac{1}{4}s^4  J_L^2 k^2 (\omega_0 \cdot L_0) \Delta{t}^4 \\
&= s^2 c^2 (J_k - J_L) k^2 \Delta{t}^2 + \frac{1}{4} s^4 \left(\Omega^2 L_0^2 - k^2 \right) J_L^3 L_0^2 \Delta{t}^4 \\
&= s^2 c^2 (J_k - J_L) k^2 \Delta{t}^2 + \frac{1}{4} s^4 \left(\Omega \cdot L_0\right)^2 J_L^3 L_0^2 \Delta{t}^4 \\
&= s^2 c^2 (J_k - J_L) k^2 \Delta{t}^2 
\end{align*}
$$
Whew. That's what we were shooting for though! This method is good because the quadratic term is *not* always positive, and as an added bonus the higher-order terms vanish (apart of course from the terms arising due to the sinc and cosine terms, which is a multiplier of $O(1 - \Omega^2 \Delta{t}^2 / 24)$) 

The next question is of course: Can we do any better? Can we nix the quadratic term entirely?

### Finding a better $J'$

We've now established that stepping the Euler equations using a rotation around a modified angular momentum $\Omega = \omega_0 - J' L_0$ can work well, showing that $J' = J_L$ produces results comparable to the Jolt method. But could a different value of $J'$ produce even better results?

There's a few different approaches we could take to this, and they all give the same answer. Let's think about the time derivatives of the rotational energy. Rather than Euler's equations, I'll treat the evolution as a rotation about $\Omega$ and we'll see what happens
$$
\begin{align*}
0 = 2 \frac{dE}{dt} &= \frac{d}{dt} (L J L) \\
&= L J \frac{dL}{dt} + \frac{dL}{dt} J L \\
&= 2 \omega \cdot \frac{dL}{dt} \\
&= 2 \omega \cdot \Omega\times L \\
&= -2 \Omega \cdot k
\end{align*}
$$
This is trivially satisfied given the form we've chosen for $\Omega$: $(\omega_0 - J' L_0)\cdot k = 0 - 0 J'$. Good start.

The second derivative gives us:
$$
\begin{align*}
2 \frac{d^2 E}{dt^2} &= L J \frac{d^2L}{dt^2} + 2 \frac{dL}{dt} J \frac{dL}{dt} + \frac{d^2L}{dt^2} J L \\
&= 2 L J \frac{d^2L}{dt^2} + 2 \frac{dL}{dt} J \frac{dL}{dt} \\
&= 2 \omega \cdot [\Omega \times (\Omega \times L)] + 2 [\Omega \times L] J [\Omega \times L] \\
&= 2 \omega \cdot (\Omega \times k) + 2 k J k \\
&= 2 [k \cdot (\omega \times \Omega) + k J k) \\
&= 2 [k J k - J' k^2]
\end{align*}
$$
There we go: For this second derivative to vanish, we need $J' k^2 = kJk$, so $J' = J_k$.

Looking back at the energy conservation expressions from earlier, we now have:
$$
\begin{align*}
2\Delta{E} &= s^2 (c^2 J_k - J') k^2 \Delta{t}^2 - s^3 c (\Omega \cdot L) (\omega_0 J k) \Delta{t}^3 + \frac{1}{4}s^4 (\Omega \times k)J(\Omega \times k) \Delta{t}^4 \\
&= s^2 (c^2 - 1) J_k k^2 \Delta{t}^2 + O(\Delta{t}^3) \\
&= -s^2 \sin^2\left(\frac{\Omega \Delta{t}}{2}\right) J_k k^2 \Delta{t}^2 + O(\Delta{t}^3) \\
&= -\frac{1}{4} s^4 J_k k^2 \Omega^2 \Delta{t}^4 + O(\Delta{t}^3) \\
&= - s^3 c (\Omega \cdot L) (\omega_0 J k) \Delta{t}^3 + \frac{1}{4}s^4 \left[(\Omega \times k)J(\Omega \times k) - J_k k^2 \Omega^2\right] \Delta{t}^4 \\
&= - s^3 c (\Omega \cdot L) (\omega_0 J k) \Delta{t}^3 + \frac{1}{4}s^4 \left[(\Omega \times k)J(\Omega \times k) - J_k \left[|\Omega \times k|^2 + (k\cdot  \Omega)^2\right]\right] \Delta{t}^4 \\
&= - s^3 c (\Omega \cdot L) (\omega_0 J k) \Delta{t}^3 + \frac{1}{4}s^4 (\Omega \times k)[J-J_k](\Omega \times k) \Delta{t}^4
\end{align*}
$$
That looks promising, we're now third-order. 

As a further note, I'll calculate here what $J_k$ looks like for a symmetric top, with eigenvalues $J_\perp, J_\perp, J_\parallel$:
$$
\begin{align*}
\omega_0 &= \left[J_\perp L_x, J_\perp L_y, J_\parallel L_z\right] \\
k &= \left[(J_\perp L_y L_z - J_\parallel L_z L_y), (J_\parallel L_z L_x - J_\perp L_x L_z), (J_\perp L_x L_y - J_\perp L_y L_x)\right] \\
&= (J_\parallel - J_\perp) L_z \left[-L_y, L_x, 0\right] \\
k^2 &= (J_\parallel - J_\perp)^2 L_z^2 \left[L_y^2 + L_x^2\right] \\
kJk &= (J_\parallel - J_\perp)^2 L_z^2 J_\perp \left[L_y^2 + L_x^2\right] \\
J_k &= J_\perp
\end{align*}
$$
Notice that, under my proposed scheme, $R(t+\Delta{t}) = Q(J' L \Delta{t}) Q(\Omega \Delta{t}) R(t)$, the finding that $J_k = J_\perp$ means we exactly reproduce the analytic integration of the symmetric top! Looking back at the energy equation to confirm for this case:
$$
\begin{align*}
J - J_k &= \mathrm{diag} [0,0,J_\parallel - J_\perp] \\
\omega_0 J k &\propto [J_\perp^2 L_x, J_\perp^2 L_y, J_\parallel L_z] \cdot [-L_y, L_x, 0] = J_\perp^2 (L_y L_x - L_x L_y) = 0 \\
\Omega \times k &\propto [0, 0, J_\parallel L_z] \times [-L_y, L_x, 0] = [-J_\parallel L_x L_z, -J_\parallel L_z L_y, 0] \\
[J - J_k] (\Omega \times k) &= 0 \\
2\Delta{E} &= - s^3 c (\Omega \cdot L) (\omega_0 J k) \Delta{t}^3 + \frac{1}{4}s^4 (\Omega \times k)[J-J_k](\Omega \times k) \Delta{t}^4 \\
&= 0
\end{align*}
$$
As expected, exact conservation for the symmetric top!

### Going Further

At this point I'm going to step beyond anything I'd necessarily suggest anyone actually implement, and talk about how we can improve on the above result. One thing we can do is write a condition relating $\Omega$ to both $L(t)$ and $L(t + \Delta{t})$. Then the hope is that we can expand that condition as a power series in $\Delta{t}$ and get higher-order $J'(\Delta{t})$. Let's try it:

One thing that we know about rotating $L(t)$ about some scaled axis $\Omega \Delta{t}$ to find $L(t + \Delta{t})$ is that the dot product between a vector and the axis of rotation is invariant under the rotation. That is:
$$
L(t) \cdot \Omega = L(t + \Delta{t}) \cdot \Omega
$$
Rearrange this to get at $J'$:
$$
\begin{align*}
[L(t + \Delta{t}) - L(t)] &\cdot (\omega_0 - J' L_0) = 0 \\
J' &= \frac{\omega_0 \cdot \Delta{L}}{L_0 \cdot \Delta{L}} 
\end{align*}
$$
Now, one thing we can do is expand $\Delta{L}$ as a power series in time. Unlike above, where we expanded the effect of rotating $L_0$ around $\Omega$, we now use Euler's equations to generate derivatives (because the effect of the rotation is already encoded in the condition that the dot product remains unchanged). In both cases we're trying to get the solution to Euler's equations and the rotation around $\Omega$ to coincide at the time $t + \Delta{t}$, but whereas before we expressed an invariant of Euler's equations ($E$) as a series in evolution by rotation, now we express an invariant of the rotation ($L \cdot \Omega$) as a series in evolution by Euler's equations. I start from the standard Taylor expansion formula and compute some derivatives (all evaluated at $t$):
$$
\begin{align*}
\Delta{L} &= \frac{dL}{dt} \Delta{t} + \frac{1}{2} \frac{d^2L}{dt^2} \Delta{t}^2 + \frac{1}{6} \frac{d^3L}{dt^3} \Delta{t}^3 + O(\Delta{t}^4) \\
\frac{dL}{dt} &= - \omega_0 \times L_0 = -k \\
\frac{d^2L}{dt^2} &= \omega_0 \times k + (Jk) \times L_0 = [(\omega_0 \cdot L_0) \omega_0 - \omega_0^2 L_0] + (Jk) \times L_0 \\
\end{align*}
$$
Now, plug the Taylor series into our expression for $J'$. To quadratic order in time, we get a constant $J'$:
$$
\begin{align*}
J' &\approx \frac{\omega_0 \cdot (-k) + \frac{1}{2} \omega_0 \cdot [(\omega_0 \cdot L_0) \omega_0 - \omega_0^2 L_0 + (Jk) \times L_0]}{L_0 \cdot (-k) + \frac{1}{2} L_0 \cdot [(\omega_0 \cdot L_0) \omega_0 - \omega_0^2 L_0 + (Jk) \times L_0]} \\
&= \frac{\omega_0 \cdot (Jk) \times L_0}{[(\omega_0 \cdot L_0)^2 - \omega_0^2 L_0^2 ]} = \frac{-kJk}{-k^2} = J_k\\
\end{align*}
$$
Well what do you know, turns out it's $J_k$ again! Let's try going one step further now:
$$
\begin{align*}
-\frac{d^3L}{dt^3} &= \frac{d^2\omega}{dt^2} \times L_0 + 2 \frac{d\omega_0}{dt} \times \frac{dL_0}{dt} + \omega_0 \times \frac{d^2L_0}{dt^2} \\
&= J\left[[(\omega_0 \cdot L_0) \omega_0 - \omega_0^2 L_0] + (Jk) \times L_0\right] \times L_0 \\
&\quad + 2 (-Jk) \times (-k) \\
&\quad + \omega_0 \times \left[[(\omega_0 \cdot L_0) \omega_0 - \omega_0^2 L_0] + (Jk) \times L_0\right] \\
&= - L_0 \times \left[[(\omega_0 \cdot L_0) J\omega_0 - \omega_0^2 \omega_0] + J[(Jk) \times L_0]\right] \\
&\quad + 2 Jk \times k \\
&\quad + \left[-\omega_0^2 k + \omega_0 \times[(Jk) \times L_0]\right] \\
&= - \left[[(\omega_0 \cdot L_0) L_0 \times J\omega_0 + \omega_0^2 k] + L_0 \times J[(Jk) \times L_0]\right] \\
&\quad - 2 (\omega_0 Jk) L_0 \\
&\quad + \left[-\omega_0^2 k + (\omega_0\cdot L_0) (Jk) - (\omega_0 Jk) L_0\right] \\
&= - \left[[(\omega_0 \cdot L_0) L_0 \times J\omega_0 ] + L_0 \times J[(Jk) \times L_0]\right] \\
&\quad - 3 (\omega_0 Jk) L_0 \\
&\quad + \left[-2\omega_0^2 k + (\omega_0\cdot L_0) (Jk)\right] \\
\omega_0 \cdot \frac{d^3L}{dt^3} &=-(\omega_0\cdot L_0) \omega_0 J k - [Jk \times L_0] \cdot Jk - 3 (\omega_0 Jk) (\omega_0 \cdot L_0) - 2 \omega_0^2 (\omega_0 \cdot k) + (\omega_0\cdot L_0) \omega_0 J k \\
&= - 3 (\omega_0 \cdot L) (\omega_0 J k) \\
L_0 \cdot \frac{d^3L}{dt^3} &= - 3 (\omega_0 J k) L_0^2 
\end{align*}
$$
Now,
$$
\begin{align*}
J' &\approx \frac{\omega_0\cdot\frac{d^2L}{dt^2} + \frac{1}{3} \omega_0\cdot\frac{d^3L}{dt^3} \Delta{t}}{L_0\cdot\frac{d^2L}{dt^2} + \frac{1}{3} L_0\cdot\frac{d^3L}{dt^3} \Delta{t}} \\
&= \frac{kJk + (\omega_0\cdot L_0)(\omega_0Jk)\Delta{t}}{k^2 + L_0^2(\omega_0Jk)\Delta{t}} \\
&= \frac{J_k k^2 + J_L L_0^2 (\omega_0Jk)\Delta{t}}{k^2 + L_0^2(\omega_0Jk)\Delta{t}} \\
&= J_k + [J_L - J_k] \frac{\alpha \Delta{t}}{1 + \alpha \Delta{t}}
\end{align*}
$$
The next-lowest order approximation to $J'$ starts at $J_k$ and decays towards $J_L$ with initial rate $\alpha = (\omega_0 J k) L_0^2 / k^2$. Makes sense that $\omega_0 Jk$ would feature, as it's part of the cubic term in the energy equation above.

You could continue this process with higher and higher derivatives of $L$, but the third is already pretty gnarly, so I won't brave that here. The one other thing I will return to briefly is the analytic solution that I mentioned near the start of this article. The analytic solution expresses the individual components $L_x, L_y, L_z$ in body space, using the Jacobi elliptic functions $\mathrm{cn}$, $\mathrm{sn}$ and $\mathrm{dn}$. These functions have addition formulae analogous to (though more complicated than) those for $\sin$ and $\cos$. This means it should be possible to express $L_{(x,y,z)}(t+\Delta{t})$ in terms of $L_{(x,y,z)}(t)$, $J_{(x,y,z)}$  and $\mathrm{cn}(\omega_p \Delta{t})$, $\mathrm{sn}(\omega_p \Delta{t})$, $\mathrm{dn}(\omega_p \Delta{t})$ (where $\omega_p$ is some given frequency dependent on the constants of motion). That would allow writing out 
$$
J' = \frac{\omega_0 \cdot (L(t + \Delta{t}) - L(t))}{L_0 \cdot (L(t + \Delta{t}) - L(t))}
$$
explicitly, which may or may not end up being a neat formula in the end. 

### My recommendation

I have two recommendations for authors of rigid body simulations, depending on their preference.

#### "Modified Jolt Method"

Continue to use the typical split between (angular) velocity integration and position (orientation) integration: 
$$
\begin{align*}
\omega_1 &= f(\omega_0, R_0) \\
R_1 &= Q(\omega_1 \Delta{t}) R_0
\end{align*}
$$

For the velocity update, use:
$$
\begin{align*}
L_0 &= J^{-1} \omega_0 + \tau \Delta{t} \\
k &= \omega_0 \times L_0 \\
J' &= \frac{k\cdot Jk}{k^2} \\
\Omega &= \omega_0 - J' L_0 \\
L_1 &= Q(-\Omega \Delta{t}) L_0 \\
\omega_1 &= J L_1
\end{align*}
$$

Optionally, use the "fast length-preserving rotation" found at the start of this post to implement the angular momentum integration:
$$
\begin{align*}
\vec{\phi} &= -\frac{1}{2} \vec{\Omega}\Delta{t} \\
u &= \phi \times L_0 \\
L_1 &= L_0 + 2 \frac{u + \phi \times u}{1 + |\phi|^2}
\end{align*}
$$


#### "Momentum-first" method

Treat the velocity update stage as a *momentum* update:
$$
\begin{align*}
L' &= L_0 + \tau \Delta{t} \\
\omega' &= (R_0 J R_0^T) L_1
\end{align*}
$$
For the position update, perform the relevant rotation and update the angular velocity to reflect the new orientation afterwards.
$$
\begin{align*}
k &= \omega' \times L' \\
J' &= \frac{k \cdot Jk }{k^2} \\
\Omega &= \omega' - J' L' \\
R_1 &= Q(J' L' \Delta{t}) Q(\Omega \Delta{t}) R_0 \\
\omega_1 &= (R_1 J R_1^T) L'
\end{align*}
$$
This treats angular momentum and orientation as the "fundamental" quantities, and angular velocity as something derived from a given momentum and orientation. Unlike the first option, it naturally preserves the direction of the angular momentum, not just its magnitude. It is a bit more of a radical departure from typical design of a rigid body integrator though.
