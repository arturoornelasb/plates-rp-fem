# E17 continuation -- 2048-window certification gates (RESULTS)

Executed 2026-07-09 on the cached production eigenpairs (no new
production solves). Two INDEPENDENT gates per sector:

| sector | two-mesh (102,63), 1.08x (original) | two-mesh (92,57), 1.196x (honest) | SSSS-exact anchor @ (110,68) | N=1024 window needs | N=2048 window needs |
|---|---|---|---|---|---|
| ee | 951 | 766 | 774 | 614 | 1229 |
| eo | 986 | 680 | 742 | 614 | 1229 |
| oe | 983 | 742 | 752 | 614 | 1229 |
| oo | 952 | 684 | 777 | 614 | 1229 |

Anchor max relative eigenvalue errors: 1.0-1.1e-3 at mode 1408.

## Reading

- The honest separated two-mesh check and the absolute exact anchor
  AGREE (680-766 vs 742-777): production statistics-grade coverage is
  ~700-770 modes/sector. The original 1.08x check (951-986) was
  correlated-optimistic, as suspected -- a documented methods point
  (close meshes share error and overstate N*).
- **The N = 1024 window ([409, 614] per sector) is certified by BOTH
  gates in ALL four sectors, with >= 66 modes of margin** -- the E17q
  verdict (PROTOCOL ARTIFACT through the certified range) stands under
  the strictest available gates.
- The N = 2048 window (needs 1229/sector) is NOT certifiable at this
  production mesh; the 2048 row remains INDICATIVE (it is flat and
  identical to the certified levels). Full 2048 certification would need
  a finer production mesh (~(140, 86) quarter, ~190k dofs/sector) --
  registered, low priority given the uniform flatness.
