
<h5>Distance to the top atmosphere boundary</h5>

<p>A point at distance $d$ from $\bp$ along $[\bp,\bi)$ has coordinates
$[d\sqrt{1-\mu^2}, r+d\mu]^\top$, whose squared norm is $d^2+2r\mu d+r^2$.
Thus, by definition of $\bi$, we have
$\Vert\bp\bi\Vert^2+2r\mu\Vert\bp\bi\Vert+r^2=r_{\mathrm{top}}^2$,
from which we deduce the length $\Vert\bp\bi\Vert$:
*/

Length DistanceToTopAtmosphereBoundary(
    Length r, Number mu) {
  assert(r <= top_radius);
  assert(mu >= -1.0 && mu <= 1.0);
  Area discriminant = r * r * (mu * mu - 1.0) +
      top_radius * top_radius;
  return ClampDistance(-r * mu + SafeSqrt(discriminant));
}

/*
<p>We will also need, in the other sections, the distance to the bottom
atmosphere boundary, which can be computed in a similar way (this code assumes
that $[\bp,\bi)$ intersects the ground):
*/

Length DistanceToBottomAtmosphereBoundary(
    Length r, Number mu) {
  assert(r >= bottom_radius);
  assert(mu >= -1.0 && mu <= 1.0);
  Area discriminant = r * r * (mu * mu - 1.0) +
      bottom_radius * bottom_radius;
  return ClampDistance(-r * mu - SafeSqrt(discriminant));
}
/*
<h4 id="transmittance_precomputation">Precomputation</h4>

<p>The above function is quite costly to evaluate, and a lot of evaluations are
needed to compute single and multiple scattering. Fortunately this function
depends on only two parameters and is quite smooth, so we can precompute it in a
small 2D texture to optimize its evaluation.

<p>For this we need a mapping between the function parameters $(r,\mu)$ and the
texture coordinates $(u,v)$, and vice-versa, because these parameters do not
have the same units and range of values. And even if it was the case, storing a
function $f$ from the $[0,1]$ interval in a texture of size $n$ would sample the
function at $0.5/n$, $1.5/n$, ... $(n-0.5)/n$, because texture samples are at
the center of texels. Therefore, this texture would only give us extrapolated
function values at the domain boundaries ($0$ and $1$). To avoid this we need
to store $f(0)$ at the center of texel 0 and $f(1)$ at the center of texel
$n-1$. This can be done with the following mapping from values $x$ in $[0,1]$ to
texture coordinates $u$ in $[0.5/n,1-0.5/n]$ - and its inverse:
*/

Number GetTextureCoordFromUnitRange(Number x, int texture_size) {
  return 0.5 / Number(texture_size) + x * (1.0 - 1.0 / Number(texture_size));
}

Number GetUnitRangeFromTextureCoord(Number u, int texture_size) {
  return (u - 0.5 / Number(texture_size)) / (1.0 - 1.0 / Number(texture_size));
}

/*
<p>Using these functions, we can now define a mapping between $(r,\mu)$ and the
texture coordinates $(u,v)$, and its inverse, which avoid any extrapolation
during texture lookups. In the <a href=
"http://evasion.inrialpes.fr/~Eric.Bruneton/PrecomputedAtmosphericScattering2.zip"
>original implementation</a> this mapping was using some ad-hoc constants chosen
for the Earth atmosphere case. Here we use a generic mapping, working for any
 but still providing an increased sampling rate near the horizon.
Our improved mapping is based on the parameterization described in our
<a href="https://hal.inria.fr/inria-00288758/en">paper</a> for the 4D textures:
we use the same mapping for $r$, and a slightly improved mapping for $\mu$
(considering only the case where the view ray does not intersect the ground).
More precisely, we map $\mu$ to a value $x_{\mu}$ between 0 and 1 by considering
the distance $d$ to the top atmosphere boundary, compared to its minimum and
maximum values $d_{\mathrm{min}}=r_{\mathrm{top}}-r$ and
$d_{\mathrm{max}}=\rho+H$ (cf. the notations from the
<a href="https://hal.inria.fr/inria-00288758/en">paper</a> and the figure
below):

<svg width="505px" height="195px">
  <style type="text/css"><![CDATA[
    circle { fill: #000000; stroke: none; }
    path { fill: none; stroke: #000000; }
    text { font-size: 16px; font-style: normal; font-family: Sans; }
    .vector { font-weight: bold; }
  ]]></style>
  <path d="m 5,85 a 520,520 0 0 1 372,105"/>
  <path d="m 5,5 a 600,600 0 0 1 490,185"/>
  <path d="m 60,0 0,190"/>
  <path d="m 60,65 180,-35"/>
  <path d="m 55,5 5,-5 5,5"/>
  <path d="m 55,60 5,5 5,-5"/>
  <path d="m 55,70 5,-5 5,5"/>
  <path d="m 60,40 a 25,25 0 0 1 25,20" style="stroke-dasharray:4,2;"/>
  <path d="m 60,65 415,105"/>
  <circle cx="60" cy="65" r="2.5"/>
  <circle cx="240" cy="30" r="2.5"/>
  <circle cx="180" cy="95" r="2.5"/>
  <circle cx="475" cy="170" r="2.5"/>
  <text x="20" y="40">d<tspan style="font-size:10px" dy="2">min</tspan></text>
  <text x="35" y="70" class="vector">p</text>
  <text x="35" y="125">r</text>
  <text x="75" y="40">μ=cos(θ)</text>
  <text x="120" y="75">ρ</text>
  <text x="155" y="60">d</text>
  <text x="315" y="125">H</text>
</svg>

<p>With these definitions, the mapping from $(r,\mu)$ to the texture coordinates
$(u,v)$ can be implemented as follows:
*/

vec2 GetTransmittanceTextureUvFromRMu(
    Length r, Number mu) {
  assert(r >= bottom_radius && r <= top_radius);
  assert(mu >= -1.0 && mu <= 1.0);
  // Distance to top atmosphere boundary for a horizontal ray at ground level.
  Length H = sqrt(top_radius * top_radius -
      bottom_radius * bottom_radius);
  // Distance to the horizon.
  Length rho =
      SafeSqrt(r * r - bottom_radius * bottom_radius);
  // Distance to the top atmosphere boundary for the ray (r,mu), and its minimum
  // and maximum values over all mu - obtained for (r,1) and (r,mu_horizon).
  Length d = DistanceToTopAtmosphereBoundary( r, mu);
  Length d_min = top_radius - r;
  Length d_max = rho + H;
  Number x_mu = (d - d_min) / (d_max - d_min);
  Number x_r = rho / H;
  return vec2(GetTextureCoordFromUnitRange(x_mu, TRANSMITTANCE_TEXTURE_WIDTH),
              GetTextureCoordFromUnitRange(x_r, TRANSMITTANCE_TEXTURE_HEIGHT));
}
