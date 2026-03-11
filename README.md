# LRTV-PDS

Low-Rank + Total Variation minimization solved via Primal–Dual Splitting (PDS) with step size adaptation.

---

## Overview

LRTV-PDS solves a convex optimization problem that combines:

- Low-rank regularization (global structure)
- Total Variation (TV) regularization (local smoothness and edge preservation)
- Noise inequality constraint
- Box constraint (range constraint on pixel values)

This framework is effective for image and tensor completion and denoising.

---

## Optimization Problem (with Box Constraint)

We solve the following convex optimization problem:

$$
\begin{aligned}
\min_{X} \quad & \alpha R_{\mathrm{rank}}(X)+ (1-\alpha) \mathrm{TV}(X) \\\\
\text{s.t.} \quad & \lVert P_{\Omega}(X - T) \lVert_p^p \le \varepsilon, \\\\
& v_{\min} \leq X \leq v_{\max} .
\end{aligned}
$$

where

- $R_{\mathrm{rank}}(X)$ is a low-rank promoting regularizer
  (e.g., tensor nuclear norm),
- $\mathrm{TV}(X)$ is total variation,
- $\alpha \in [0,1]$ can be tuned manually,
- $P_{\Omega}$ is the observation operator (for removing missing entries),
- $T$ is the observed tensor,
- $\varepsilon$ controls the noise level,
- $p=1$ or $p=2$ can be selected,
- $v_{\min}$ and $v_{\max}$ are minimum and maximum values for all pixels.

---


## Applications

- Image inpainting  
- Image denoising  
- Tensor completion  
- Video restoration  

---

## Citation

If you use this code, please cite:

```bibtex
@article{yokota2019simultaneous,
  title={Simultaneous tensor completion and denoising by noise inequality constrained convex optimization},
  author={Yokota, Tatsuya and Hontani, Hidekata},
  journal={IEEE Access},
  volume={7},
  pages={15669--15682},
  year={2019},
  publisher={IEEE}
}
```
