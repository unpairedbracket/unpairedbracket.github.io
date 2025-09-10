---
tags: 
  - physics
  - game design
  - character control
---

<style>
    iframe {
        width: 70%; 
        height: auto;
        aspect-ratio: 16/9;
        margin: auto;
        display: block;
    }
</style>

> [!note] Work-in-Progress
> This post is under construction, read at your own peril!

I've been thinking about designing a character controller for a game project. I want nebulously "physics-based" platforming, with a high-speed character whose speed is affected by things like running up (slower) or down (faster) slopes. In addition, when moving sufficiently fast running on vertical surfaces, or ceilings should be possible, *and control intuitively*.

## Basic motion

On the most basic level, let's start with kinematics "on the flat" and in one dimension. I want my character to have a maximum speed, and a maximum acceleration. Then there's a sort of "characteristic time" to the motion, given by the ratio of these quantities:
$$
t_\mathrm{acc} = \frac{v_\mathrm{max}}{a_\mathrm{max}}
$$
Consider motion to the right ($+x$), at full throttle (analog stick fully to the right, or just right key pressed on a keyboard). Acceleration should start at $a_\mathrm{max}$ if I'm not moving and reduce to zero if I reach $v_\mathrm{max}$.
$$
a(v) = \frac{\partial v}{\partial t} = \frac{v_\mathrm{max} - v}{t_\mathrm{acc}}
$$
This equation of motion will produce an exponential velocity profile, as a function of time:
$$
v(t) = v_\mathrm{max} \left(1 - \exp(-t / t_\mathrm{acc})\right)
$$
This simple formula can also be modelled as applying a constant acceleration force equal to $m v_\mathrm{max} / t_\mathrm{acc}$ against a velocity-proportional drag force $-m v / t_\mathrm{acc}$. The mass of the character will be largely ignored in this kinematic treatment, but keep in mind that in the context of a real game, it would affect the forces the character imparts on other objects, e.g. when it runs into things.

<iframe src="char_motion_example/side_scroll_basic.html"></iframe>

### Slowing down
Once we're moving at speed, we might at some point want to slow down or stop. Either we release the "go forwards" control, or switch to a "go back" one. 

With no control direction applied, we could just trust the drag force from above to slow us down. I find this feels "slippery" though, so apply a braking enhancement, which increases the drag force by reducing the effective acceleration time:
$$
t_\mathrm{stop} = t_\mathrm{acc} / \alpha_\mathrm{stop}
$$
$$
\frac{\partial v}{\partial t} = - \frac{v}{t_\mathrm{stop}}
$$
I've found a $\alpha_\mathrm{stop} = 5$ to be a good value, leading to our character slowing down five times faster than it accelerates up to speed. 

<iframe src="char_motion_example/side_scroll_stopping.html"></iframe>

### Reversing direction
What if we're moving to the right and start pushing left? When the control direction and the current velocity are in opposite directions, we enhance the "drag" component of the force by a factor $\alpha_\mathrm{rev}$
$$
\frac{\partial v}{\partial t} = \frac{v_\mathrm{targ} - \alpha_\mathrm{rev} v}{t_\mathrm{acc}} 
$$
This can be thought of as modifying the acceleration time, only when the control direction and current velocity are in opposition.
$$
t_\mathrm{rev} = t_\mathrm{acc} \frac{v_\mathrm{targ} - v}{v_\mathrm{targ} - \alpha_\mathrm{rev} v}
$$
Note also that this looks similar to the slowing to a stop case, but with an applied target velocity and potentially a different slowing factor

<iframe src="char_motion_example/side_scroll_reverse.html"></iframe>

### Going faster than top speed
At some points, through interaction with the game world or objects, e.g. going down slopes, the character may end up moving faster than its maximum speed. I want to allow the player to go fast, and reward them by allowing them to continue going fast once they're there. The current model will result in decay of velocity back towards the maximum with the base acceleration time scale $t_\mathrm{acc}$. To allow the player to maintain high speeds for longer, a "slowing factor" $\beta_\mathrm{slow}$ is applied when the character is moving faster than the maximum speed - but of course, only when the control direction is parallel to the current velocity.

## The flat, one-dimensional model

This completes the model for one dimensional motion:
$$
v_\mathrm{targ} = c v_\mathrm{max}
$$
$$
\frac{\partial v}{\partial t} = \frac{v_\mathrm{targ} - \alpha v}{\tau}
$$

$$
\alpha = \begin{cases}
\alpha_\mathrm{stop} & c = 0 \\
\alpha_\mathrm{rev} & v\cdot c < 0 \\
1 & \mathrm{otherwise}
\end{cases}
$$
$$
\tau = \begin{cases} 
\beta_\mathrm{slow} t_\mathrm{acc} & v\cdot \hat{c} > v_\mathrm{max} \\ 
t_\mathrm{acc} & \mathrm{otherwise}
\end{cases}
$$
This model trivially generalises to analogue input using an input intent variable $-1 < c < 1$: now $v_\mathrm{targ}$ includes a factor of $c$ so can take any value between $-v_\mathrm{max}$ and $+v_\mathrm{max}$. $\hat{c}$ is the input direction, i.e. in one dimension either $-1$, $+1$ or zero.

## A second dimension, taken one way: top-down motion

The simplest way to add in a second dimension is to expand our ground from a flat line to a flat plane. Now the control direction goes from being a single number between $-1$ and $1$ to a unit (or zero) vector, and our character can move freely within the plane.

The main difference from the one-dimensional case here is that we can now _turn_. Start with the basic force-drag equation of motion established before:
$$
\frac{\partial \vec{v}}{\partial t} = \frac{\vec{c}\ v_\mathrm{max} - \vec{v}}{t_\mathrm{acc}}
$$
This works fine, but as with stopping and reversing direction in one dimension, it feels "slippery" on its own. Turning is as gradual as acceleration, so changing direction feels very sluggish. The solution we used earlier for stopping is good for the velocity component transverse to the control direction:
$$
\begin{align*}
\frac{\partial v_\parallel}{\partial t} &= \frac{v_\mathrm{targ} - v_\parallel}{t_\mathrm{acc}} \\
\frac{\partial v_\perp}{\partial t} &= -  \frac{\alpha_\mathrm{turn} v_\perp}{t_\mathrm{acc}}
\end{align*}
$$
We've aligned coordinates according to the control direction, such that $(c v_\mathrm{max})_\parallel = |c| v_\mathrm{max} = v_\mathrm{targ}$ and $(c v_\mathrm{max})_\perp = 0$

There's a problem here: Because transverse velocity decays to zero faster than longitudinal velocity evolves, sudden changes in direction end up looking weird. If the component of velocity along the control axis ($v_\parallel$) is negative, the stronger transverse acceleration can cause it to "snap" to moving in the direction *opposite* to the control direction, then more gradually slow to a stop and start accelerating towards the control direction.

<iframe src="char_motion_example/top_down_problem.html"></iframe>

To solve this problem, I modify the longitudinal acceleration when $v_\parallel$ is negative. The acceleration should be in between one directly towards the target velocity, and one towards a complete stop. At any given moment, the point along the $v_\parallel$ axis that the acceleration points towards is:
$$
v' = v_\parallel - \frac{a_\parallel}{a_\perp} v_\perp
$$
To fulfil $0 < v' < v_\mathrm{targ}$ impose the condition:
$$
\begin{align*}
0 &< v_\parallel + \frac{a_\parallel t_\mathrm{acc}}{\alpha_\mathrm{turn}} < v_\mathrm{targ} \\
0 &< \alpha_\mathrm{turn} v_\parallel + a_\parallel t_\mathrm{acc} < \alpha_\mathrm{turn} v_\mathrm{targ}
\end{align*}
$$
I also want to match the forward acceleration value at $v_\parallel = 0$:
$$
t_\mathrm{acc} a_\parallel(v_\parallel = 0) = v_\mathrm{targ} - v_\parallel = v_\mathrm{targ}
$$
There's a space of solutions here, but the one I choose is: $a_\parallel t_\mathrm{acc} = v_\mathrm{targ} - \alpha_\mathrm{turn} v_\parallel$. This matches the approach to changing direction I showed in the one-dimensional case above (as long as we unify the $\alpha_\mathrm{turn}$ and $\alpha_\mathrm{rev}$ coefficients), which is nice (and no coincidence)!

When stopping ($\vec{c} = 0$), the assignment of parallel and transverse components breaks down. In that case, I apply the stopping enhancement with a factor $\alpha_\mathrm{stop}$ to the entire velocity vector, matching the behaviour in one dimension. I've generally found that using the same value ($5$, in particular) for all three $\alpha$ parameters produces a nice snappy feeling to character motion, but your mileage may vary of course. Play around!

<iframe src="char_motion_example/top_down.html"></iframe>

## A second dimension, the other way: Side-scrolling

Now I'll change perspective, to a different two dimensions. I'll retain one horizontal direction for the character to move along, but the second dimension will now be "up". The new, interesting thing we can now do is allow the ground to slope. Gently - not so much that the character is unable to stand still or walk in either direction; slipping and sliding will be handled later.

When the character is moving downhill its maximum speed should be increased, and when moving uphill it should be decreased. Because we still want to be able to move in either direction, the maximum speed should not be decreased to zero or become negative. 

The ground is sloped by an angle $\theta$, between $-\pi/2$ (vertically down) and $\pi/2$ (vertically up). I don't allow overhangs yet, that's a subject for a different blog post. Naïvely, an acceleration (along the sloped ground) of $g_\mathrm{slope} = -g\sin{\theta}$ would be added to the equation of motion. However, gravity stronger than $a_\mathrm{max}$ could then remove the ability of the player to walk uphill or even to stand still.

I'd like for steep slopes to be able to more than double the character's maximum speed while maintaining the ability to move uphill, albeit slowly. To achieve this I choose to pass the "physical" effect of gravity through a filter, to reduce the effect of moving uphill and increase the effect of moving downhill. Rather than the physically-based modification to maximum speed:
$$
v_\mathrm{slope} = v_\mathrm{flat} + g_\mathrm{slope} t_\mathrm{acc} = v_\mathrm{flat} \left(1 + \frac{g_\mathrm{slope}}{a_\mathrm{max}} \right)
$$
I use a "softplus"-style function (this particular functional form will be better-motivated later on, when I move back to top-down motion and introduce 3D):
$$
v_\mathrm{slope} = v_\mathrm{flat} \left(\sqrt{1 + \gamma^2 \frac{g_\mathrm{slope}}{a_\mathrm{max}}^2} + \gamma \frac{g_\mathrm{slope}}{a_\mathrm{max}}\right)
$$
with a tweakable parameter $\gamma$ to allow adjusting the strength of the slope effect independently of gravity itself (within what I've discussed in this post, there's no point in this, but it might be helpful for tuning the slope effect independently of, say, jump behaviour). $\gamma = 1$ replicates the more physical behaviour for small slopes, $\gamma = 0.5$ keeps the effect of downhill gravity mostly similar to the physical model for large values of $g_\mathrm{slope}$. To make the effect of gravity feel more dramatic I have found values between $\gamma=1$ and $\gamma=10$ to work well.

<iframe src="char_motion_example/side_scroll_hills.html"></iframe>

## Synthesis: All Three Dimensions

Just as we've added a vertical dimension to the initial one-dimensional model to arrive at a side-scrolling 2D model, I'll now consider adding a vertical dimension to the top-down 2D controller described above. The challenge is similar: slopes and how they affect the motion[^control_plane]. As before, I want the following things:
- Character can stand still even on arbitrarily steep slopes, with arbitrarily strong gravity
- The steeper a slope, the lower the maximum speed when running uphill and the faster you can go downhill
- The strength of the slope effect should be greater downhill than uphill - make it easier to go fast, and keep going fast!
I also want:
- The maximum component of the character's speed along a contour line - the direction orthogonal to up/downhill to be their "normal" top speed, the same as on the flat
  - By this I don't mean that motion along a contour (i.e. on a slope, but "across" the slope rather than upwards or downwards) should be the same speed as on the flat, but that when running downhill at some angle to the slope, regardless of how fast you're going down the slope your maximum component of speed across the slope should be the normal max velocity.
- The *direction* of the target velocity shouldn't be altered by the slope, only its magnitude. If I push along a contour line I don't want to drift downhill due to gravity.

So I need a function that scales the character's maximum velocity as a function of the downhill gravity component and the angle between the control direction and the downhill direction.

It won't be immediately apparent that the solution presented here is a generalisation of the "side scroller" solution I showed above, but I'll show how it reduces at the end.

To start with, now that our ground is a plane, $\vec{g}_\mathrm{slope}$ becomes a vector - the projection of the gravitational acceleration onto the ground plane. Call the ground plane normal $\hat{n}$; then $\vec{g}_\mathrm{slope} = \vec{g} - (\hat{n} \cdot \vec{g})\hat{n}$.

The function I'll use will be an ellipse, with one focus at the origin of "velocity scaling space" and the other at
$$
2\frac{\beta\vec{g}_\mathrm{slope}}{a_\mathrm{max}}
$$
Its semi-minor axis will be 1, to maintain a constant maximum along-the-contour velocity component. 

The equation of an ellipse, in polar coordinates relative to one focus, is:
$$
r(\theta) = \frac{b \sqrt{1-e^2}}{1 - e \cos{\theta}}
$$
where $e$ is the eccentricity, $b$ the semi-minor axis and $\theta = 0$ in the direction toward the other focus: $\hat{c} \cdot \vec{g}_\mathrm{slope} = |g_\mathrm{slope}| \cos\theta$.

The linear eccentricity $c$ of an ellipse is the distance between its centre and either focus, and is therefore
$$\frac{\beta\vec{g}_\mathrm{slope}}{a_\mathrm{max}}$$

From standard ellipse relations:
$$
\begin{align}
c = ae \\
c^2 = a^2 - b^2
\end{align}
$$
with $a$ the semi-major axis. Eliminating $a$ and solving for $e$,
$$
\begin{align}
c^2 &= c^2 / e^2 - b^2 \\
e^2 (b^2 + c^2) &= c^2 \\ 
e &= \frac{c}{\sqrt{b^2 + c^2}}
\end{align}
$$

This gives a "vector eccentricity" of: 
$$
\vec{e} = \frac{\beta\vec{g}_\mathrm{slope}}{\sqrt{a_\mathrm{max}^2 + \beta^2g_\mathrm{slope}^2}}
$$
which has the nice property that the $e \cos\theta$ in the denominator of the ellipse equation is expressible as $\vec{e} \cdot \vec{c}$.

The numerator of the ellipse equation is the "semi-latus rectum"[^no_laughing], notable as the speed multiplied purely along a contour line (at $\theta = \pi / 2$). This is less than $v_\mathrm{max}$, which is the largest the contour line component can be when moving at any angle to the slope - it turns out it's maximised at $\cos\theta = |e|$.

This speed factor at $\theta = \pi / 2$ is:
$$
\ell = b \sqrt{1 - e^2} = \frac{a_\mathrm{max}}{\sqrt{a_\mathrm{max}^2+\beta^2g_\mathrm{slope}^2}}
$$
And so the overall factor is:
$$
r(\hat{c}) = \frac{a_\mathrm{max}}{\sqrt{a_\mathrm{max}^2+\beta^2g_\mathrm{slope}^2} - \beta \hat{c} \cdot \vec{g}_\mathrm{slope}}
$$

When travelling perfectly downhill ($\hat{c} \parallel +\vec{g}_\mathrm{slope}$) or uphill ($\hat{c} \parallel -\vec{g}_\mathrm{slope}$) this reduces to:
$$
r^\parallel= \frac{a_\mathrm{max}}{\sqrt{a_\mathrm{max}^2+\beta^2g_\mathrm{slope}^2} \mp \beta |g|_\mathrm{slope}} = \frac{\sqrt{a_\mathrm{max}^2+\beta^2g_\mathrm{slope}^2} \pm \beta |g|_\mathrm{slope}}{a_\mathrm{max}} = \sqrt{1+\left(\frac{\beta |g|_\mathrm{slope}}{a_\mathrm{max}}\right)^2} \pm \frac{\beta |g|_\mathrm{slope}}{a_\mathrm{max}}
$$
which is, as promised, equivalent to the "side-scrolling" slope handling result presented above

[^control_plane]: I also need to define how to map the control direction onto the ground plane, but that's a story for another day. 
[^no_laughing]: Stop laughing.