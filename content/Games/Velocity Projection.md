We're designing a move-and-slide system for a physics engine.
The basic idea is to shapecast the body in a given direction, by a given distance (i.e. $v \Delta t$). If an object is hit by the shapecast, we stop at the hit location (after some time $\Delta t_\mathrm{hit}$), adjust the velocity to point in a direction that wouldn't move inside the hit object, and continue with a second shapecast, using the new velocity vector and the remaining time. This process is iterated until a shapecast completes without hitting anything.

The question is how to do the "adjust the velocity to point in a direction that wouldn't move into the hit object" bit. For a single colliding object, this is pretty easy and obvious: look at the collision normal $\hat{n}$ (defined to be the direction _away_ from the colliding plane, the "hat" indicates that it's a unit vector) and the velocity $\vec{v}$ (the arrow indicating a vector that's not necessarily a unit vector). Then the "not moving into the hit object" condition is expressed as $\hat{n} \cdot \vec{v} \geq 0$.

If this is already fulfilled (if you're standing right next to a wall and start moving away from it this could happen) that's fine, return the velocity as it is. If it's not, the simplest thing is just to zero out the velocity component pointing into the hit plane: $\vec{v}_\mathrm{next} \cdot \hat{n} = 0$ which we can achieve by saying $\vec{v}_\mathrm{next} = \vec{v} - (\vec{v}\cdot\hat{n}) \hat{n}$. Because $\hat{n}$ is a unit vector, $\vec{v}_\mathrm{next}\cdot\hat{n} = \vec{v}\cdot\hat{n} - \vec{v}\cdot\hat{n} = 0$.

Okay, that wasn't too hard. What if we're in contact with more than one thing? Say there are now two planes we can't move into. We now have more than one condition that must be fulfilled simultaneously: $\hat{n}_j\cdot\vec{v} \geq 0$ for $j$ in $\{1,2\}$. I'll accumulate the normals into a matrix, $N = \left[\hat{n}_1, \hat{n}_2\right]$ so I can write this as $N^T \vec{v} \geq \vec{0}$, where the inequality between vectors is defined elementwise. I'll call these equations, one per plane normal, _the constraints_.

As before, if the input velocity is consistent with both of the plane constraints, nothing needs to be done and we can return it unaltered. If only one of the constraints is violated, it looks like the process used for a single-plane projection can be used on that single constraint and will return a valid result. This case is shown in the first diagram below. 

![[2_normals.svg]]If both constraints are violated, we could try solving each individual constraint in turn. However, as shown in the middle diagram, the result is not independent of the order in which we evaluate the constraints. Evaluating the constraint arising from normal pointing down and left first produces the yellow result, which then fulfils the other constraint without needing further changes. If we evaluate the "down and right" constraint first, the "down and left" constraint remains violated and solving that constraint then leaves us in a different place to the yellow result. This sequence of orthogonal projections is shown in orange. Maybe, if (as in the yellow case) resolving any one constraint resolves all of the others incidentally, we should prioritise that one?

There isn't always a single constraint that can resolve all others though: what should be done in the third diagram? Resolving the two constraints sequentially is again order-dependent, and this time there isn't an order that results in only having to resolve a single constraint.

It's also possible to come up with pairs of constraints for which this sequential resolution process doesn't result in a valid output _at all_. What should we do there? Resolving such a situation using sequential projections might require repeated alternating projections, repeating until some degree of convergence is found:
![[converging.svg]]
It seems like we need some quality metric to let us say which is the "best" point on the feasible set. A reasonable one is that the impulse imparted on the moving body by the hit objects is as small in magnitude as possible, that is: $\min |\vec{v}_\mathrm{out} - \vec{v}_{in}|^2$. Then we can say that the yellow solution in the second diagram, and the cyan one in the third, look like the best feasible solutions.

We now have a proper problem statement:

$$\vec{v}_\mathrm{out} = \min_{\vec{v}} |\vec{v} - \vec{v}_{in}|^2 \quad \mathrm{s.t.}\quad N^T \vec{v} \geq 0$$

This is a quadratic programming problem, and amounts to orthogonal projection onto the "cone" of feasible directions defined by the constraints. This "cone" is, in the two-dimensional drawings above, the unshaded wedge towards the bottom of the circle
![[voronoi.svg]]
Here, for each of the two cases shown above, "Voronoi regions" for the full range of velocity directions are shown. The green region shows the directions of velocities that will not be changed by the constraints: they already fulfil all of them. The closest point in the unshaded wedge to any of the green points is the green point itself. The red regions towards the tops of the circles are closest to the point at the centre, i.e. zero velocity. Notice that in each case the red region is the region between the points where the two normals would intersect the circle if traced backwards from the origin - these "negative normals" are shown in dotted lines. The blue regions of the circle are closest to one of the radii that define the allowed sector, velocities in those directions should be projected onto the corresponding straight edge of the unshaded wedge. You can see from this diagram that one way of solving the problem is to calculate the shortest distance between the input velocity and:
 - Each of the sides of the unshaded wedge
	 - This shortest distance is given by the dot product of the constraint normal and the input velocity (blue regions)
 - The origin 
	 - In this case the distance is the magnitude of the input velocity (red region)
We also need to check if each of the nearest points is actually valid according to all of the constraints, and check for the case where all constraints are already met by the input (green region).
## Three dimensions
In three dimensions (unfortunately it'll be harder to draw things now), three planes can form three intersections. Imagine looking at the corner of a cube, from the inside - if you're indoors look up at the top corner of the room you're in. You'll see three planes (the ceiling and two walls) which correspond to the blue lines in the diagrams from the previous section. The red normals are the normals to the walls/ceilings that point into the room, and the blue shaded regions are on the sides of the walls outside of the room you're in. There are three edges, one between each pair of the two walls and the ceiling, and they all meet at the corner. As you're inside the room, you're in the feasible region of the constraints already - a three-dimensional wedge of space corresponding to the unshaded sector in the previous section. This means an object in contact with the corner would be able to move towards you unrestricted.
If such an object were trying to move towards someone in the next room over, on the other side of one of the walls, the closest feasible point would be their closest point on the wall between them and you. The object would slide along the dividing wall.
Similarly, the closest feasible point for someone standing in the room above you is the point they're standing on, so to try and move towards them the object would slide along the ceiling.
If someone is in the room above one of those neighbouring rooms, _their_ closest feasible point is somewhere along the seam between their floor (your ceiling) and the wall they're on the other side of, and someone on the floor above and on the other side of both walls has theirs in the corner you're looking at itself. Trying to move toward the first of these people would result in sliding along in contact with both the ceiling and the wall; movement toward the second is impossible.

The key difference between this and the 2D case is the presence of intersections between the constraint planes, which gives us three new things to check against.

So for a given input velocity direction, we now need to check the shortest distance from that point in velocity space to:
- The three normal constraints (walls/ceilings)
- The three intersections between them (seams between walls/ceilings)
- The corner itself
We take whichever of these is closest, and that's our solution. We still also need to check that all of our constraints are satisfied at each point, and reject it if they're not.

## Many Constraints

In some cases, with complex static colliders for example, there may be many constraints to solve. The number of intersections between $N$ constraint planes is $N(N-1)/2$. The algorithm described above now has to deal with:
 - The $N$ individual constraints
 - The $N(N-1)/2$ intersections between constraint planes
 - The "corner" at the origin
Including the $N-1$ checks to check whether a projection onto a constraint plane is valid, and the $N-2$ necessary per intersection, and $N$ initial checks to determine whether the initial velocity is legal, we're now having to perform potentially  O(N^3)$ tests. That's not ideal.

So what can we do? For a particular set of $N$ constraint normals, some number $M \leq N$ actually define the polyhedral cone of allowed velocities (corresponding to the unshaded region on the 2d diagrams above), and those $M$ limiting constraint normals are exactly the ones that fall on the exterior of the convex hull of normal vectors (defined on $\mathbb{S}^2$ in the sense of geodesic convexity[^1], not as a conventional convex hull in $\mathbb{R}^3$) - this corresponds to the red Voronoi region shown above.

More helpfully, only $M$ of the total $N(N-1)/2$ intersections between constraint planes can ever be the limiting point for an input velocity, and these correspond to pairs of normals that share an edge in the aforementioned convex hull. The other intersections are always forbidden by unrelated constraints.
In other words, adding a new constraining normal to a system has no effect, _if_ that new normal falls within the convex hull of the existing set of normals

The new difficulty is that it is not straightforward to know - in general - which loop of $M$ normals from a given set of $N$ form the convex hull. Indeed it is not even trivial to determine the value $M$ itself. In spherical geometry and under the geodesic definition of convexity, we must also handle the complication that the convex hull of a set of points only _has_ a boundary if those points all fall on some given hemisphere. Otherwise, the convex hull is the entire sphere. Determining whether this is the case is also not trivial. All of this means it's quite difficult to determine which constraint (or pair of constraints) are the important ones to resolve.

## Dual Problem
Let's look at the problem from a different perspective. Recall the quadratic programming problem we're trying to solve (negating the constraint condition to put the result in what is apparently a more "standard" form):

$$\vec{v}_\mathrm{out} = \arg\min_{\vec{v}} \frac{1}{2}|\vec{v} - \vec{v}_\mathrm{in}|^2 \quad \mathrm{s.t.}\quad -N^T \vec{v} \leq 0$$

I'm going to introduce a vector $\vec{\lambda}$ of $N$ non-negative Lagrange multipliers and derive the [Wolfe dual](https://en.wikipedia.org/wiki/Wolfe_duality) (the Lagrangian dual problem works out the same but this is simpler to work through, or at least makes more sense to me):
$$
\begin{align}
\max_{\vec{\lambda}\geq 0,\vec{v}} \frac{1}{2}&|\vec{v} - \vec{v}_\mathrm{in}|^2 + \vec{\lambda} \cdot(-N^T\vec{v}) \\
\mathrm{s.t.}\quad \nabla_{v}\frac{1}{2}&|\vec{v} - \vec{v}_\mathrm{in}|^2 - \nabla_{v} (N\vec{\lambda}) \cdot \vec{v} = 0
\end{align}
$$
Evaluating the gradients in the constraint gives
$$
\begin{align}
(\vec{v} &- \vec{v}_\mathrm{in}) - N\vec{\lambda} = 0 \\
\vec{v} &= \vec{v}_\mathrm{in} + N\vec{\lambda}.
\end{align}
$$
This allows us to make intuitive sense of the Lagrange multipliers: each one is the strength of the "impulse" (velocity change) imparted by its corresponding constraining plane. We can define the total velocity impulse: $\Delta{\vec{v}} = N\vec{\lambda}$.
Substitute the expression for $\vec{v}$ into the function we're trying to maximise:
$$
\begin{align}
\vec{\lambda} = \arg\max_{\vec{\lambda}\geq 0} &\frac{1}{2}|N\vec{\lambda}|^2 - \vec{v} \cdot N \vec{\lambda} \\
= &\frac{1}{2}N\vec{\lambda} \cdot (N\vec{\lambda} - 2(\vec{v}_\mathrm{in} + N\vec{\lambda})) \\ 
= &- \frac{1}{2}N\vec{\lambda} \cdot (2\vec{v}_\mathrm{in} + N\vec{\lambda})) \\ 
= &- N\vec{\lambda} \cdot \frac{\vec{v}_\mathrm{in} + \vec{v}}{2} \\ 
= &- (\vec{v} - \vec{v}_\mathrm{in}) \cdot \frac{\vec{v}_\mathrm{in} + \vec{v}}{2} \\
= & \frac{|v_\mathrm{in}|^2 - |v|^2}{2}
\end{align}
$$
### Constraining the output velocity
As an aside, because this is a fun fact about the problem: Because quadratic optimisation with linear constraints exhibits _strong duality_, we know that - at the optimum point - the values of the primal problem's objective function and that of this new dual problem are equal.

Therefore,

$$\frac{1}{2} |\vec{v} - \vec{v}_\mathrm{in}|^2 = \frac{1}{2} |N\vec{\lambda}|^2 - \vec{v} \cdot N\vec{\lambda}.$$

From the constraint equation we evaluated just now, we know that $\vec{v} - \vec{v}_\mathrm{in} = N\vec{\lambda}$, so the first term on the right-hand side is the same as the entire left-hand side. Cancelling gives $\vec{v} \cdot N\vec{\lambda} = 0$, or in other words $\vec{v} \cdot \Delta{\vec{v}} = 0$. The final velocity is orthogonal to the change in velocity! Because $\vec{v} = \vec{v}_\mathrm{in} + \Delta{\vec{v}}$, this orthogonality implies that the final velocity and the change in velocity form a right-angled triangle whose hypotenuse is the input velocity!
Now, to finally, for the first time in my life, draw on a circle theorem I learned in secondary school many years ago, and subsequently had no use for over the course of a cumulative 10 years as an undergraduate and then PhD student in physics[^2]: If you draw a line across the diameter of a circle, and then connect both ends to any point on that circle, the resulting lines form a right-angle.
Looking at it the other way around: all possible right-angled triangles with a given line segment as their hypotenuse have their right angle touching the circle (or sphere in 3d) of which the hypotenuse is a diameter. Because of the way spheres work, I can also express the result as $\vec{v} = \frac{1}{2}\vec{v}_\mathrm{in} + \frac{1}{2}|v_\mathrm{in}| \hat{a}$. Then, $N\vec{\lambda} = \frac{1}{2}(|v_\mathrm{in}| \hat{a} - \vec{v}_\mathrm{in})$. I don't think this helps us at all, but I didn't try very hard to make use of it.
![[circle_theorem.svg]]
The upshot of this is that we can draw a diagram like this: initial velocity in green and as before, normals in red and forbidden regions shaded in blue. The difference is we're now drawing the normal cone as extending from $\vec{v}_\mathrm{in}$ rather than from the origin. The solution is the point on the perimeter of the circle that's closest to the top (to $\vec{v}_\mathrm{in}$) and not in a blue region, or equivalently the lowest point (closest to $\vec{0}$) that _is_ between the red normals. It's obvious that the two agree.
### back to actually solving this thing
I'll now return to the dual problem I derived above, before I got distracted. Dropping constant offsets and scalings from the expression to be maximised, and flipping the sign to change from a maximisation to a minimisation, we have:

$$\vec{\lambda} = \arg\min_\mathrm{\lambda} |\vec{v}_\mathrm{in} + N\vec{\lambda}|^2 \quad \mathrm{s.t.}\quad \vec{\lambda} \geq 0$$

In plain English: We're looking for a combination of the constraint normals, with non-negative coefficients $\vec{\lambda}$, that minimises the magnitude of $\vec{v}_\mathrm{in} + N \vec{\lambda}$. Writing out the origin $O$ explicitly, $|(\vec{v}_\mathrm{in} + N \vec{\lambda}) - \vec{0}| = d(O, \vec{v}_\mathrm{in} + N \vec{\lambda})$, the distance between the output velocity and the origin. Rearranging though, it's also $|\vec{0} + N\vec{\lambda} - (-\vec{v}_\mathrm{in}) | = d(O + N \vec{\lambda}, -\vec{v}_\mathrm{in})$. An optimal set of coefficients $\vec{\lambda}$ represents the point of the _normal cone_ (NC) that's closest to the negative of the initial velocity. The velocity-space coordinates of that point on the NNC are $N\vec{\lambda} = \Delta{\vec{v}}$ , so the output velocity, $\vec{v} = \Delta{\vec{v}} - (-\vec{v}_\mathrm{in})$ is the vector pointing to that closest point from the negative of the input velocity.

This is a point projection onto convex polytope problem, exactly the kind that's solved by the Gilbert-Johnson-Keerthi algorithm. Unfortunately, the convex polytope is an infinite polyhedral cone, which I would not expect to be well-supported by any existing GJK implementation. 

## GJK on $\mathbb{S}^2$

My initial idea was to try and run GJK on the unit sphere. Place the origin at the centre of the unit sphere and rotate space to point the initial velocity upwards, then use a Gnomonic projection to map the normal cone down onto a plane. This projection would project the negative input velocity to the origin of the plane, map the negative normals to points and the geodesic arcs between them to straight lines. It doesn't preserve spherical distances _in general_, but distances _from the origin_ of the plane are a monotonic function of the corresponding spherical distances. All this means that if the closest point to the origin is a vertex of the convex hull in the plane, then the closest point of the normal cone on the sphere is the corresponding normal, and if the closest point in the plane is on the line between two of the vertices, then the closest point of the normal cone is on the arc between the corresponding normals. Once you know which normal or two define the projected point, projecting the velocity onto those one or two constraints is straightforward. Finally, If the negative input velocity is found _within_ the convex hull of the NC, then it trivially _is_ the closest point within the NC to itself.

The problem here is that the Gnomonic projection is does not work for vectors pointing in arbitrary directions. If a normal is pointing *away* from the projection plane ($\hat{n}\cdot\vec{v}_\mathrm{in} > 0$), projecting it does not work in a way that's consistent with those projecting _towards_ it. Ultimately I was not able to find a straightforward and efficient way to handle such normals, which cut short my attempts to solve the problem using weird sphere-plane projections (at least, after thinking about stereographic projection for a bit and deciding it wasn't going to help me).

## GJK Against the Conical Hull

At this point I considered returning to three dimensions and trying to clip the infinite cone to finite extents, but sufficiently large to not cause problems. This didn't seem like it would be possible though, because depending on the arrangement of the cone normals the required lengths could be arbitrarily large. The lengths needed are also dependent on _every other_ normal in the constraint set, adding a $O(N^2)$ factor in that I'd prefer to do without. 

So my next thought was: What if I do just let the cone be infinite (at least, arbitrarily large)? The support function of such a shape is actually pretty simple: the support point in any given direction is represented by the normal whose dot product with that point is greatest, unless they're all negative in which case the support point is the apex of the cone (the origin).

It's easy to seed the process with the nearest point on the cone to the input: the origin is the only  point for which the distance to an input is finite. The closest edge of the cone is also simple (and $O(N)$) to find: the support normal in the direction of the input velocity (from the origin) defines the ray that passes closest to the input velocity.

The infinite analogs to the simplices of GJK are:
 - Point -> The Origin (this is the only point not infinitely far from the input point)
 - Edge -> Ray (between the origin and a point at infinity)
 - Triangle -> Wedge (the area between two rays, the virtual third edge of the triangle is infinitely far from the input point)
 - Tetrahedron -> Solid Wedge (bounded by three Rays and the three Wedges between them; the fourth face of the tetrahedron is at infinity)

Because the origin and the initial ray are the closest point and edge respectively to the target point, subsequent iterations need only check for proximity of simplex faces (Wedges) to the target. This means that for a triangle it is only necessary to find the projected point on the face (one vertex is the origin, the other two are at infinity; two sides are rays from the origin to infinity and the other is wholly at infinity. If the origin or the closer of the two rays contained the goal point, the algorithm would have terminated already, leaving only the face to evaluate). For a tetrahedron it is only necessary to check the two _new_ faces (the other one was checked in the previous stage, and the arguments about rays and points translate directly from the triangle case).
This enables some optimisations in the "walking tetrahedra" stage of GJK.

The overall algorithm goes as follows:
 - Input the list of normals $\hat{n}_j$ and the initial velocity, and set $\vec{x}_0 = -\vec{v}_\mathrm{in}$.
 - Find the index $j$ for which $\vec{x}_0 \cdot \hat{n}_j$ is maximised. The ray formed by the origin and the normal $\hat{n}_j$ is the closest ray to $x_0$. The closest point to $\vec{x}_0$ along the ray is $(\vec{x}_0 \cdot \hat{n}_j) \hat{n}_j$.
 - The "search direction" is the vector $\vec{s}_1 = \vec{x}_0 - (\vec{x}_0 \cdot \hat{n}_j)\hat{n}_j$ from the closest point on ray $j$ towards $\vec{x}_0$. So: find the normal $n_k$ with the maximum dot product $\vec{s}_1 \cdot \hat{n}_k$. Add this normal to the "simplicial cone", which is now the wedge between the origin and rays $j$ and $k$. The closest point on this wedge must be on its surface. If it were on ray $j$, we would have terminated on the previous step. If it were on ray $k$, that ray would have been the first one selected as it would be closer than ray $j$. The two rays span a plane whose normal is proportional to $\vec{e}_{jk} = \hat{n}_j \times \hat{n}_k$. The closest point to $\vec{x}_0$ is now $\vec{x}_0 - (\hat{e}_{jk} \cdot \vec{x}_0) \hat{e}_{jk}$, where $\hat{e}_{jk}$ is the normalisation of $\vec{e}_{jk}$. That means the new search direction is $\vec{s}_2 = (\hat{e}_{jk} \cdot \vec{x}_0) \hat{e}_{jk}$. 
 - Once again, find the normal with the greatest dot product against $\vec{s}_2$ and call it $\hat{n}_l$. Again, if no normal with positive dot product exists, return the current best point. Add this to our wedge to form a... Solid wedge? Not sure what's best to call it. It's like a tetrahedron with one vertex at the origin and the opposite face blown out to infinity. One of the wedge faces is the one we've already analysed, so the closest point isn't on that one. It can't be on the new ray, for the reason described above. Now we've got a solid wedge, $\vec{x}_0$ could be inside the volume of the wedge without being on its boundary. For each of the two new wedge faces, check if $\vec{x}_0$ is on the same side of the wedge as the opposing normal is (i.e. is the sign of $n_j \cdot \vec{e}_{kl}$ the same as the sign of $\vec{x}_0 \cdot \vec{e}_{kl}$ and so on) - we already know this is the case for the old wedge and the new normal, so just need to check the other two. If that is the case for both new wedges as well, then $\vec{x}_0$ is inside the current solid wedge and the algorithm can terminate. Otherwise, retain the normals that form the wedge with the greatest signed distance to $\vec{x}_0$ (oriented such that the distance from the wedge to its opposing normal is considered negative) and repeat this step with this new wedge and new search direction $\vec{s}_\mathrm{new} = (\hat{e}_\mathrm{new}\cdot\vec{x}_0)\hat{e}_\mathrm{new}$.
 - The above step repeats, "walking" along solid wedges until one a) contains $\vec{x}_0$ or b) cannot be expanded in the direction of $\vec{x}_0$, meaning the global best point is on its boundary. 
 - The found "closest point on normal cone to $\vec{x}_0$" is $N\vec{\lambda}$. This needs to be added to $\vec{v}_\mathrm{in}$ to yield the output velocity. $\vec{x}_\mathrm{best} + \vec{v}_\mathrm{in} = \vec{x}_\mathrm{best} - \vec{x}_0$ is the negation of the final search direction (which is zero if the search terminated by finding $\vec{x}_0$ within the normal cone).

Having implemented and tested this strange version of GJK it seems to be both fast and robust. Although the runtime is in principle $O(N^2)$ for $N$ normals (as each step involves an $O(N)$ search for the next supporting normal, and up to $N$ steps may be taken in total), I have found it difficult to construct sets of normals for which the algorithm takes more than four steps to terminate for any given input velocity.

Here's an attempt:
![[5_steps_clipped.svg]]
From left to right:
 - The negated input velocity is the green point at the centre of the circle. The closest normal to the input is the red point. No normal can ever be found closer than that, marked by the large red shaded circle. The search direction is to the left, shown in blue
 - A second normal is found to the left. No other normal can now be found further to the left; the forbidden region is marked with blue shading. The new search direction is, again, towards the input. Again, shown with the blue arrow.
 - A third normal now forms a triangle on the surface. Once again, all points further along the previous search direction than the new one are now forbidden. The new search direction is in blue, and the allowed region is now quite small
 - The triangle is flipped to include a new, fourth normal within the small remaining allowed region. Another region is forbidden, and the remaining allowed region (circled in green) is so small it's covered by the strokes of the lines.
 - Zooming in, the thin sliver of space in which to put a fifth normal becomes visible. By continuing this procedure, you can imagine adding normals exponentially-smaller regions of space. Practically though, it looks like more than four iterations will be exceedingly rare.
 If you find a reliable construction which leads to more than four steps for a significant area of the sphere, do let me know!

## Returning to the primal problem

You may be wondering how this relates to the original problem we were trying to solve. Basically it turns out that the "support point" I defined for the polyhedral cone corresponds to finding the _most strongly violated constraint_ in the primal problem. So thinking about the problem in the primal description again, the algorithm looks something like this:
 - For a given input velocity, find the most strongly violated constraint (most negative value of $\hat{n}_j \cdot \vec{v}$). If all of them are non-negative then no constraints are violated and the input velocity is accepted immediately. Otherwise, resolve this constraint by orthogonal projection onto the constraint plane.
 - From the _new_ velocity, after resolving the first constraint, find the new most strongly violated one. Again, if no constraints are now violated we can accept the new velocity as a solution. Otherwise, project the input velocity onto the intersection of the first active constraint and the new one
 - With this new velocity, find the new most strongly violated constraint. If there aren't any, terminate and accept the current projection. If there are, add that one to the set. You now have three active constraints: one new one and two old ones. We know that simply resolving the new one on its own won't work, because if it did it would have been selected in the first step and the algorithm would have terminated. Calculate the projection of the input velocity onto the intersections between the new constraint plane and each of the old ones, and see if each projection is allowed according to the one currently active constraint it's _not_ projected onto. If both are forbidden, we can kill all velocity because there's no valid point other than zero velocity. If one is allowed and the other forbidden, keep the two normals that define the allowed intersection and drop the remaining one. If _both_ projections are allowed, go with the one that results in projected velocity of greatest magnitude. Iterate this step until all constraints are satisfied, or until the "kill all velocity" condition is met.
I find it quite hard to translate thinking in terms of the dual problem to think in terms of the primal problem again, so some of the details above might be wrong. I think the basic gist and structure are right though. In any case, I think implementing a solution for the dual problem is the neater and easier to visualise option.

[^1]: The spherical convex hull can be thought of as covering every point on the sphere you can normalising a weighted sum (in $\mathbb{R}^3$) of the input points, where the weights are non-negative. That weighted sum is known as a _conical combination_ of the inputs. The condition seen in Euclidean convex hulls that the weights sum to 1 is replaced by a condition that the end result has a _length_ of 1, which is achieved by the normalisation. 

[^2]: Not sure if this is a win or a loss for the people who say "you'll never need to know this outside of the exam" about the maths we learn in school. Closer to a win than usually, I guess
