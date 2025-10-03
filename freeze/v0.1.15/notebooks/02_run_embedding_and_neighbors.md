# 02 — Run embedding & neighbors (scripted)

This step is scripted. Use the Makefile targets; no notebook logic lives here.

```bash
# Build/refresh processed objects & models (idempotent)
make prep

# Attach .raw & (optional) symbol mapping
make map

# Compute neighbors/UMAP on X_scvi
make umap
