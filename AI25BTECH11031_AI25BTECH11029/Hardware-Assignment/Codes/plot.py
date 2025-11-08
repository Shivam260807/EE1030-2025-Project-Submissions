import numpy as np
import matplotlib.pyplot as plt

# 3D figure setup
fig = plt.figure(figsize=(8, 6))
ax = fig.add_subplot(111, projection='3d')

# Weighted inner product matrix (diagonal)
D = np.diag([1, 0.5, 1/3])

# Vector for polynomial (1 + x)
u = np.array([1, 1, 0])

# Effective normal in Euclidean space (weighted inner product)
n = D @ u   # n = (1, 0.5, 0)

# Plot the vector (1+x)
ax.quiver(0, 0, 0, u[0], u[1], u[2],
          color='red', linewidth=2, label='Vector $(1+x)$')

# Define the orthogonal plane: n·x = 0 → a0 + 0.5a1 = 0 → a0 = -0.5a1
a1 = np.linspace(-2, 2, 20)
a2 = np.linspace(-2, 2, 20)
A1, A2 = np.meshgrid(a1, a2)
A0 = -0.5 * A1  # plane equation

# Plot the orthogonal plane
ax.plot_surface(A0, A1, A2, color='lightblue', alpha=0.5)

# Basis vectors of the orthogonal subspace
v1 = np.array([-0.5, 1, 0])  # corresponds to (1 - 2x)
v2 = np.array([0, 0, 1])     # corresponds to (x^2)

# Plot basis vectors
ax.quiver(0, 0, 0, v1[0], v1[1], v1[2],
          color='green', linewidth=2, label='Basis $(1-2x)$')
ax.quiver(0, 0, 0, v2[0], v2[1], v2[2],
          color='purple', linewidth=2, label='Basis $(x^2)$')

# Labels and style
ax.set_xlabel('$a_0$')
ax.set_ylabel('$a_1$')
ax.set_zlabel('$a_2$')
ax.set_title('Orthogonal Subspace to $(1+x)$ under Weighted Inner Product')
ax.legend()
ax.view_init(elev=25, azim=45)

plt.show()

