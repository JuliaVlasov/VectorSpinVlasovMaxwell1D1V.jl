export Advection

Advection(mesh::Mesh, t::BSpline) = BSplineAdvection(mesh, t.p)
Advection(mesh::Mesh, t::PSM) = PSMAdvection(mesh)
