#include "cfd/rf/network.hpp"
#include "cfd/rf/nport.hpp"
#include "cfd/rf/antenna.hpp"
#include "cfd/rf/peec.hpp"
#include "cfd/rf/thin_wire_mom.hpp"
#include "cfd/circuit/spice.hpp"
#include "cfd/circuit/analysis.hpp"
#include "cfd/multiphysics/block_system.hpp"
#include "cfd/multiphysics/deformation_optics.hpp"
#include "cfd/fvm/combustion.hpp"
#include "cfd/fvm/temporal.hpp"
#include "cfd/fvm/radiation_network.hpp"
#include "cfd/solvers/fvm/premixed_flame1d.hpp"
#include "cfd/fvm/function_objects.hpp"
#include "cfd/solvers/fdtd/cylindrical_tm.hpp"
#include "cfd/core/backend.hpp"
#include "cfd/core/profiler.hpp"
#include "cfd/core/block_csr.hpp"
#include "cfd/core/task_graph.hpp"
#include "cfd/core/simd.hpp"
#include "cfd/core/scratch_arena.hpp"
#include "cfd/core/mixed_precision.hpp"
#include "cfd/core/geometric_multigrid.hpp"
#include "cfd/core/numa.hpp"
#include "cfd/core/amg_adapter.hpp"
#include "cfd/core/distributed_csr.hpp"
#include "cfd/fvm/reordering.hpp"
#include "cfd/solvers/lbm/visualization.hpp"
#include "cfd/core/csr_matrix.hpp"
#include "cfd/fvm/field.hpp"
#include "cfd/fvm/runtime_schemes.hpp"
#include "cfd/fvm/poly_mesh.hpp"
#include "cfd/fem/mesh2d.hpp"
#include "cfd/multiphysics/projection.hpp"
#include "cfd/io/mesh_io.hpp"
#include "cfd/io/openfoam_io.hpp"
#include "cfd/optics/system_models.hpp"
#include "cfd/solvers/optics/materials.hpp"
#include "cfd/optics/transform.hpp"
#include "cfd/optics/analysis.hpp"
#include "cfd/optics/optimization.hpp"
#include "cfd/optics/freeform.hpp"
#include "cfd/optics/serialization.hpp"
#include "cfd/electrochemistry/corrosion_models.hpp"
#include "cfd/electrochemistry/electrochemistry.hpp"
#include "cfd/solvers/lbm/passive_scalar.hpp"
#include "cfd/solvers/fvm/thermal_transport.hpp"
#include "cfd/solvers/fvm/species_transport.hpp"
#include "cfd/solvers/fem/dynamic_bar1d.hpp"
#include "cfd/solvers/fem/brinkman2d.hpp"
#include "cfd/optics/diffraction.hpp"
#include "cfd/optics/nonsequential.hpp"
#include "cfd/workflow/schema.hpp"
#include "cfd/workflow/provenance.hpp"
#include "cfd/workflow/validation_catalog.hpp"
#include "cfd/fvm/derived_fields.hpp"
#include "cfd/solvers/lbm/diagnostics.hpp"
#include "cfd/solvers/fdtd/postprocess.hpp"
#include "cfd/chemistry/thermo.hpp"
#include "cfd/chemistry/advanced_thermo.hpp"
#include "cfd/chemistry/reactor_network.hpp"
#include "cfd/chemistry/kinetics.hpp"
#include "cfd/chemistry/equilibrium.hpp"
#include "cfd/chemistry/reaction_path.hpp"
#include "cfd/optics/wavefront.hpp"
#include "cfd/optics/birefringence.hpp"
#include "cfd/workflow/runner.hpp"
#include "cfd/fvm/fv_matrix.hpp"
#include "cfd/fvm/turbulence.hpp"
#include "cfd/solvers/fvm/vof.hpp"
#include "cfd/solvers/fvm/particles.hpp"
#include "cfd/solvers/fem/eddy_current2d.hpp"
#include "cfd/fem/material_laws.hpp"
#include "cfd/fem/phase_change.hpp"
#include "cfd/physics/material.hpp"
#include "cfd/solvers/fdtd/maxwell3d.hpp"
#include "cfd/solvers/fvm/reacting.hpp"
#include "cfd/solvers/fvm/compressible1d.hpp"
#include "cfd/fvm/advanced_models.hpp"
#include "cfd/solvers/fem/stokes2d.hpp"
#include "cfd/fem/nonlinear_geometry.hpp"
#include "cfd/physics/material_bridge.hpp"
#include "cfd/multiphysics/cht.hpp"
#include "cfd/fvm/radiation.hpp"
#include "cfd/fvm/mesh_motion.hpp"
#include "cfd/multiphysics/electrochemical_heat.hpp"
#include "cfd/solvers/lbm/advanced_d2q9.hpp"
#include "cfd/solvers/lbm/voxelize.hpp"
#include <chrono>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>
#include <numbers>
#include <stdexcept>
#include <thread>
#include <vector>
namespace {
void require(bool ok, const char* msg) {
    if (!ok)
        throw std::runtime_error(msg);
}
bool near(double a, double b, double tol = 1e-10) {
    return std::abs(a - b) <= tol * std::max({1.0, std::abs(a), std::abs(b)});
}
void test_runtime() {
    using namespace cfd::core;
    require(backend_available(RuntimeBackend::serial), "serial backend missing");
    auto b = select_backend("auto");
    require(backend_available(b), "auto backend unavailable");
    require(parse_backend("serial") == RuntimeBackend::serial, "backend parse");
    Profiler::global().reset();
    {
        ScopedProfile p("completion.unit");
        std::this_thread::sleep_for(std::chrono::milliseconds(1));
    }
    { ScopedProfile p("completion.unit"); }
    auto snap = Profiler::global().snapshot();
    require(snap.size() == 1 && snap[0].calls == 2 && snap[0].total_seconds > 0, "profiler counter");
    BlockCsrBuilder<2> builder(2, 2);
    BlockCsrMatrix<2>::Block a{2, 0, 0, 3}, b01{1, 0, 0, 1}, a11{4, 0, 0, 5};
    builder.add(0, 0, a);
    builder.add(0, 1, b01);
    builder.add(1, 1, a11);
    auto A = builder.build();
    std::vector<double> x{1, 2, 3, 4}, y(4);
    A.multiply(x, y);
    require(near(y[0], 5) && near(y[1], 10) && near(y[2], 12) && near(y[3], 20), "block CSR multiply");
}
void test_fvm_field() {
    using namespace cfd::fvm;
    DimensionedField<double> T("T", dim_temperature, FieldLocation::cell, 4, 300.0);
    require(T.size() == 4 && T.dimensions() == dim_temperature, "dimensioned field");
    FieldHistory<double> h(3);
    h.push({1, 2});
    h.push({3, 4});
    h.push({5, 6});
    require(h.old(2)[0] == 1, "field history");
    require(parse_face_interpolation_scheme("vanLeer") == FaceInterpolationScheme::van_leer, "runtime scheme");
    auto m = make_cartesian_hexa_mesh(2, 1, 1);
    std::vector<double> flux(m.face_count(), 0);
    for (std::size_t i = 0; i < m.face_count(); ++i)
        flux[i] = m.faces()[i].area.x;
    double dt = courant_limited_timestep(m, flux, 0.5, 10);
    require(dt > 0 && dt < 10, "Courant dt");
}
void test_projection() {
    using namespace cfd;
    auto fem = fem::make_rectangle_tri_mesh(8, 6, 1.0, 1.0);
    std::vector<double> nodal(fem.node_count());
    for (std::size_t i = 0; i < nodal.size(); ++i)
        nodal[i] = 2.0 + 3.0 * fem.nodes[i].x - 0.5 * fem.nodes[i].y;
    auto fvm = fvm::make_cartesian_hexa_mesh(5, 4, 1);
    auto cell = multiphysics::fem_nodes_to_fvm_cells(fem, nodal, fvm);
    for (std::size_t c = 0; c < cell.size(); ++c) {
        auto p = fvm.cells()[c].center;
        require(near(cell[c], 2 + 3 * p.x - 0.5 * p.y, 1e-9), "FEM->FVM linear projection");
    }
    std::vector<double> constant(fvm.cell_count(), 7.0);
    auto back = multiphysics::fvm_cells_to_fem_nodes(fvm, constant, fem, 6);
    for (double v : back)
        require(near(v, 7), "FVM->FEM constant projection");
    multiphysics::StructuredScalarGrid3D g;
    g.nx = g.ny = g.nz = 3;
    g.values.resize(27);
    for (std::size_t k = 0; k < 3; ++k)
        for (std::size_t j = 0; j < 3; ++j)
            for (std::size_t i = 0; i < 3; ++i)
                g.values[(k * 3 + j) * 3 + i] = double(i + j + k) / 2.0;
    require(near(g.sample(.25, .25, .25), .75), "structured trilinear sample");
}
void test_io() {
    namespace fs = std::filesystem;
    auto dir = fs::temp_directory_path() / "cfd_completion_io";
    fs::create_directories(dir);
    auto gm = dir / "tri.msh";
    {
        std::ofstream o(gm);
        o << "$MeshFormat\n2.2 0 8\n$EndMeshFormat\n$Nodes\n4\n1 0 0 0\n2 1 0 0\n3 1 1 0\n4 0 1 "
             "0\n$EndNodes\n$Elements\n2\n1 2 0 1 2 3\n2 2 0 1 3 4\n$EndElements\n";
    }
    auto m = cfd::io::read_gmsh_v2_tri2d(gm);
    require(m.node_count() == 4 && m.element_count() == 2, "Gmsh Tri3 import");
    std::vector<double> scalar{0, 1, 2, 3};
    auto vtu = dir / "mesh.vtu";
    cfd::io::write_vtu(vtu, m, scalar, "phi");
    require(fs::file_size(vtu) > 100, "VTU output");
    auto obj = dir / "mesh.obj";
    {
        std::ofstream o(obj);
        o << "v 0 0 0\nv 1 0 0\nv 0 1 0\nf 1 2 3\n";
    }
    require(cfd::io::read_obj_triangles(obj).triangles.size() == 1, "OBJ import");
    auto stl = dir / "mesh.stl";
    {
        std::ofstream o(stl);
        o << "solid x\n facet normal 0 0 1\n outer loop\n vertex 0 0 0\n vertex 1 0 0\n vertex 0 1 0\n endloop\n "
             "endfacet\nendsolid\n";
    }
    require(cfd::io::read_ascii_stl(stl).triangles.size() == 1, "STL import");
    auto round = cfd::io::read_vtu_tri2d(vtu);
    require(round.node_count() == m.node_count() && round.element_count() == m.element_count(), "VTU Tri3 roundtrip");
    auto foam = dir / "polyMesh";
    fs::create_directories(foam);
    {
        std::ofstream o(foam / "points");
        o << "8\n(\n(0 0 0)\n(1 0 0)\n(1 1 0)\n(0 1 0)\n(0 0 1)\n(1 0 1)\n(1 1 1)\n(0 1 1)\n)\n";
    }
    {
        std::ofstream o(foam / "faces");
        o << "6\n(\n4(0 3 2 1)\n4(4 5 6 7)\n4(0 1 5 4)\n4(1 2 6 5)\n4(3 7 6 2)\n4(0 4 7 3)\n)\n";
    }
    {
        std::ofstream o(foam / "owner");
        o << "6\n(\n0\n0\n0\n0\n0\n0\n)\n";
    }
    {
        std::ofstream o(foam / "neighbour");
        o << "0\n(\n)\n";
    }
    {
        std::ofstream o(foam / "boundary");
        o << "1\n(\nwalls\n{\ntype wall;\nnFaces 6;\nstartFace 0;\n}\n)\n";
    }
    auto pm = cfd::io::read_openfoam_polymesh(foam);
    require(pm.cell_count() == 1 && pm.face_count() == 6 && near(pm.cells()[0].volume, 1.0, 1e-12),
            "OpenFOAM polyMesh cube import");
    require(pm.patch_index("walls") == 0, "OpenFOAM boundary patch import");
    auto fld = dir / "T";
    cfd::io::write_openfoam_vol_scalar_field(fld, "T", {300.0});
    auto fv = cfd::io::read_openfoam_vol_scalar_field(fld, 1);
    require(fv.size() == 1 && near(fv[0], 300.0), "OpenFOAM scalar field roundtrip");
    auto uf = dir / "U";
    cfd::io::write_openfoam_vol_vector_field(uf, "U", {{1.0, 2.0, 3.0}});
    require(fs::file_size(uf) > 100, "OpenFOAM vector field output");
    fs::remove_all(dir);
}
void test_optics() {
    using namespace cfd::optics;
    auto tr = RigidTransform::from_euler_xyz({1, 2, 3}, 0.1, -0.2, 0.3);
    Ray r{{.2, .3, .4}, {0, 0, 1}, 550};
    auto rr = tr.inverse(tr.apply(r));
    require(near(rr.origin.x, r.origin.x, 1e-12) && near(rr.origin.y, r.origin.y, 1e-12) &&
                near(rr.origin.z, r.origin.z, 1e-12),
            "optical transform roundtrip");
    WavelengthModel w;
    w.add({486.1, 1, false});
    w.add({587.6, 1, true});
    require(near(w.primary_nm(), 587.6), "wavelength model");
    FieldModel f(FieldType::angle);
    f.add({0, 5, 1});
    require(f.fields().size() == 1, "field model");
    Aperture a{ApertureKind::rectangular, 2, 1, true};
    require(a.contains(1.9, .9) && !a.contains(2.1, 0), "general aperture");
    MeritFunction q = [](std::span<const double> x) { return (x[0] - 2) * (x[0] - 2) + 4 * (x[1] + 1) * (x[1] + 1); };
    auto opt = coordinate_descent(q, {0, 0}, {1, 1}, 200, 1e-8);
    require(opt.converged && std::abs(opt.variables[0] - 2) < 1e-6 && std::abs(opt.variables[1] + 1) < 1e-6,
            "local optimizer");
    auto sens = finite_difference_sensitivity(q, std::vector<double>{2.5, -.5});
    require(std::abs(sens[0] - 1) < 1e-6 && std::abs(sens[1] - 4) < 1e-6, "optical sensitivity");
    auto mc = monte_carlo_tolerance(q, std::vector<double>{2, -1}, std::vector<double>{.1, .1}, 500, 42);
    require(mc.samples == 500 && mc.mean > 0, "Monte Carlo tolerancing");
    auto de = differential_evolution(q, std::vector<double>{-5, -5}, std::vector<double>{5, 5}, 28, 120, 0.8, 0.9, 7);
    require(de.merit < 1e-6 && std::abs(de.variables[0] - 2) < 2e-3 && std::abs(de.variables[1] + 1) < 2e-3,
            "global differential evolution");
    const double derivative = autodiff_derivative([](auto x) { return square(x - decltype(x){3.0}); }, 5.0);
    require(near(derivative, 4.0, 1e-12), "optical autodiff primitive");
    MaterialCatalog catalog;
    const auto& nearest = catalog.nearest_index(587.6, 1.5168);
    require(nearest.name == "N-BK7", "glass nearest-index search");
    PolynomialFreeformSurface freeform(0.0, {{2, 0, 0.1}, {0, 2, 0.2}}, 2.0);
    Vec3 hit{}, normal{};
    Ray fr{{0.2, 0.1, -1.0}, {0, 0, 1}, 550};
    require(freeform.intersect(fr, hit, normal) && near(hit.z, 0.1 * 0.04 + 0.2 * 0.01, 1e-10),
            "freeform ray intersection");
    SequentialOpticalSystem sys(1.0);
    SequentialSurface surf;
    surf.type = SurfaceType::conic;
    surf.vertex_z = 1.0;
    surf.radius = 2.0;
    surf.aperture_radius = 1.0;
    surf.refractive_index_after = 1.5;
    surf.conic_constant = -1.0;
    sys.add_surface(surf);
    auto json = sequential_system_to_json(sys);
    auto round = sequential_system_from_json(json);
    require(round.surfaces().size() == 1 && round.surfaces()[0].type == SurfaceType::conic &&
                near(round.surfaces()[0].radius, 2.0),
            "optical JSON roundtrip");
}
void test_electrochem() {
    using namespace cfd::electrochemistry;
    RedoxActivityTerm prod[] = {{0.1, 1}}, react[] = {{1.0, 1}};
    const double E = concentration_dependent_nernst(0.5, 298.15, 1, prod, react);
    require(near(E, nernst_potential(0.5, 298.15, 1, 0.1), 1e-12), "concentration Nernst");
    auto o = oxygen_reduction_cathodic(0.2, 1.229, 0.2, 1e-7, 1.0, 20.0);
    require(o.limited_current_density <= 0 && std::abs(o.limited_current_density) <= 20.0 + 1e-12,
            "oxygen cathodic limit");
    auto h = hydrogen_evolution_cathodic(-0.2, 1.0, 1e-3, 0.5, 10.0);
    require(h.limited_current_density <= 0, "hydrogen cathodic");
    const double rct = charge_transfer_resistance(2.0);
    auto z0 = randles_impedance(0, 3, rct, 1e-3);
    require(near(z0.real(), 3 + rct, 1e-12) && near(z0.imag(), 0, 1e-12), "EIS DC");
    auto zh = randles_impedance(1e9, 3, rct, 1e-3);
    require(std::abs(zh.real() - 3) < 1e-4, "EIS high-frequency");
    DoubleLayerState dl{0.2, 0};
    dl.advance_current(2.0, 0.1);
    require(near(dl.potential, 1.0), "double layer integration");
    require(near(bruggeman_effective_transport(2.0, 0.25, 1.5), 0.25, 1e-12), "porous electrode Bruggeman transport");
    require(near(solid_phase_ohmic_drop(10.0, 0.1, 5.0), 0.2, 1e-12), "solid-phase electronic conduction");
    require(std::abs(marcus_current_density(1.0, 0.1, 50000.0)) > 0 &&
                near(marcus_current_density(1.0, 0.0, 50000.0), 0.0, 1e-12),
            "Marcus kinetics");
    require(bikerman_activity_correction(1.0, 0.1, 2.0) > 1.0, "finite-size PNP correction");
    auto iw = smooth_electrode_interface(0.0, 0.01);
    require(near(iw.electrolyte_fraction, 0.5, 1e-12) && iw.delta > 0, "immersed electrode interface");
    auto coupled = solve_anodic_bv_mass_transfer(1.0, 1e-4, 0.0, 0.08, 2.0, 298.15, 2.0, 1.0);
    require(coupled.surface_concentration > 1.0 && coupled.current_density > 0 &&
                near(1e-4 * (coupled.surface_concentration - 1.0), coupled.molar_flux, 1e-8),
            "implicit BV mass-transfer coupling");
    const std::vector<cfd::chemistry::AcidFamily> acid{{0.1, 0, {1e-5}}};
    double ep = acid_speciation_redox_potential(0.0, 298.15, 1.0, 1.0, acid, 0.05);
    require(std::isfinite(ep), "pH coupled Nernst potential");
    require(pitting_initiates(0.2, 0.5, {0.1, 0.4}), "pitting initiation criterion");
    require(stress_assisted_exchange_current(1.0, 100e6, 1e-6) > 1.0, "stress assisted corrosion");
    require(phase_field_corrosion_step(0.2, 0.0, 0.5, 1.0, 1.0, 0.01) > 0.2, "phase-field corrosion driving");
    require(near(recessed_level_set(-0.1, 0.02), -0.08, 1e-12), "moving corrosion level set");
}
void test_extended_physics() {
    using namespace cfd;
    lbm::PassiveScalarD2Q5 scalar({64, 32, 0.02, 0.0, 0.0});
    scalar.initialize([](double x, double) { return 1.0 + 0.1 * std::sin(2.0 * 3.14159265358979323846 * x); });
    const double total0 = scalar.total();
    const auto before = scalar.scalar();
    double a0 = 0;
    for (std::size_t y = 0; y < 32; ++y)
        for (std::size_t i = 0; i < 64; ++i)
            a0 += before[y * 64 + i] * std::sin(2.0 * 3.14159265358979323846 * (static_cast<double>(i) + 0.5) / 64.0);
    a0 *= 2.0 / (64.0 * 32.0);
    scalar.step(100);
    const auto after = scalar.scalar();
    double a1 = 0;
    for (std::size_t y = 0; y < 32; ++y)
        for (std::size_t i = 0; i < 64; ++i)
            a1 += after[y * 64 + i] * std::sin(2.0 * 3.14159265358979323846 * (static_cast<double>(i) + 0.5) / 64.0);
    a1 *= 2.0 / (64.0 * 32.0);
    require(std::abs(scalar.total() - total0) / total0 < 1e-12, "passive scalar conservation");
    require(a1 < a0 && a1 > 0, "passive scalar diffusion");
    auto mesh = fvm::make_cartesian_hexa_mesh(4, 4, 1);
    fvm::ThermalTransport heat(mesh, {0.01, 1.0, 2.0, 0.0});
    heat.initialize(300.0);
    heat.set_volumetric_source([](fvm::Vec3, double) { return 2.0; });
    heat.run(10);
    for (double t : heat.temperature())
        require(near(t, 300.1, 1e-11), "thermal uniform source");
    auto buoy = fvm::boussinesq_acceleration(mesh, heat.temperature(), 300.0, 0.01, {0, -9.81, 0});
    require(near(buoy[0].y, 0.00981, 1e-10), "Boussinesq source");
    fvm::MultiSpeciesTransport species(mesh, {0.01, 0.02}, 1e-3);
    species.initialize(0, [](fvm::Vec3 p) { return 1 + p.x; });
    species.initialize(1, [](fvm::Vec3 p) { return 2 - p.y; });
    double inv0 = species.inventory(0), inv1 = species.inventory(1);
    species.run(5);
    require(near(species.inventory(0), inv0, 1e-9) && near(species.inventory(1), inv1, 1e-9),
            "multi-species conservation");
    fem::DynamicBar1DConfig dc;
    dc.elements = 24;
    dc.dt = 2e-6;
    fem::DynamicBar1D bar(dc);
    const double pi = 3.14159265358979323846;
    bar.initialize([&](double x) { return 1e-4 * std::sin(0.5 * pi * x / dc.length); });
    const double e0 = bar.kinetic_energy() + bar.strain_energy();
    bar.run(100);
    const double e1 = bar.kinetic_energy() + bar.strain_energy();
    require(std::abs(e1 - e0) / e0 < 2e-4, "Newmark bar energy");
    auto bm = fem::make_rectangle_tri_mesh(10, 6);
    fem::Brinkman2D br(std::move(bm), {1.0, 1.0, 0.1, 3000, 1e-10});
    br.solve([](fem::Node2) { return fem::Displacement2{1.0, 0.0}; });
    double vmax = 0;
    for (auto v : br.velocity())
        vmax = std::max(vmax, v.x);
    require(vmax > 0, "Brinkman flow response");
    std::vector<std::complex<double>> pupil(32, {1.0, 0.0});
    auto psf = optics::normalized_psf_1d(pupil);
    double sum = 0;
    for (double v : psf)
        sum += v;
    require(near(sum, 1.0, 1e-12), "PSF normalization");
    auto mtf = optics::mtf_from_psf_1d(psf);
    require(near(mtf[0], 1.0, 1e-12), "MTF zero frequency");
    auto ff = optics::fraunhofer_1d(pupil, 1e-5, 500e-9, 1.0);
    require(ff[ff.size() / 2].intensity > 0.99, "Fraunhofer central maximum");
    optics::NonSequentialScene scene;
    scene.add_plane({{0, 0, 0}, {0, 0, -1}, 1.0, 1.0, 1.0, 0.0});
    auto paths = scene.trace({{{0, 0, -1}, {0, 0, 1}, 550}, 1.0, 0}, 2);
    require(paths.size() == 1 && paths[0].ray.direction.z < 0 && near(paths[0].power, 1.0),
            "non-sequential mirror trace");
}

void test_workflow_and_diagnostics() {
    using namespace cfd;
    workflow::CaseConfig cfg;
    cfg.set("dt", {1e-3, "s"});
    cfg.set("rho", {1000, "kg/m3"});
    cfg.validate({{"dt", "s", 1e-9, 1.0, true}, {"rho", "kg/m3", 0.0, {}, true}});
    require(cfg.to_json().find("parameters") != std::string::npos, "case schema json");
    auto parsed_json = workflow::case_config_from_json(cfg.to_json());
    require(near(parsed_json.at("dt").value, 1e-3) && parsed_json.at("dt").unit == "s", "case JSON parse");
    auto parsed_yaml = workflow::case_config_from_yaml("dt: 0.001 s\nrho: 1000 kg/m3\n");
    require(near(parsed_yaml.at("rho").value, 1000) && parsed_yaml.at("rho").unit == "kg/m3", "case YAML parse");
    workflow::MaterialDatabase db;
    db.add({"water", {{"density", {998.2, "kg/m3"}}, {"viscosity", {1e-3, "Pa.s"}}}});
    require(near(db.at("water").property.at("density").value, 998.2), "material database");
    auto prov = workflow::make_provenance("completion", "0.8.0", "test");
    prov.metadata["backend"] = "cpu";
    require(workflow::provenance_json(prov).find("git_commit") != std::string::npos, "provenance json");
    require(workflow::validation_catalog().size() >= 6, "validation catalog");
    auto mesh = fvm::make_cartesian_hexa_mesh(5, 5, 1);
    std::vector<fvm::Vec3> u(mesh.cell_count());
    for (std::size_t i = 0; i < u.size(); ++i) {
        auto p = mesh.cells()[i].center;
        u[i] = {-p.y, p.x, 0};
    }
    auto d = fvm::velocity_derived_fields(mesh, u);
    double omega = 0;
    for (auto w : d.vorticity)
        omega += w.z;
    omega /= d.vorticity.size();
    require(std::abs(omega - 2.0) < 0.25, "vorticity derived field");
    require(fvm::wall_y_plus(1e-3, 1.0, 1000, 1e-3) > 0, "y+ helper");
    std::vector<double> pressure(mesh.cell_count(), 2.0);
    auto fx = fvm::pressure_force_on_patch(mesh, pressure, mesh.patch_index("right"));
    require(fx.x < 0, "pressure patch force");
    require(near(fvm::sample_nearest_cell(mesh, pressure, {.3, .3, .3}), 2.0), "field probe");
    lbm::D2Q9Solver solver({24, 24, .7F});
    for (std::size_t y = 10; y < 14; ++y)
        for (std::size_t x = 10; x < 14; ++x)
            solver.set_solid(x, y);
    auto ft = lbm::momentum_exchange_force(solver, 12, 12);
    require(std::abs(ft.force_x) < 1e-5 && std::abs(ft.force_y) < 1e-5, "LBM equilibrium force symmetry");
    auto drag = lbm::darcy_forchheimer_acceleration(1.0, -2.0, 0.1, 0.5, 0.3);
    require(drag[0] < 0 && drag[1] > 0, "porous drag direction");
    std::vector<lbm::Vector2> vel(25);
    std::vector<double> rho(25);
    for (std::size_t y = 0; y < 5; ++y)
        for (std::size_t x = 0; x < 5; ++x) {
            vel[y * 5 + x] = {-static_cast<double>(y), static_cast<double>(x)};
            rho[y * 5 + x] = static_cast<double>(x + y);
        }
    auto vort = lbm::vorticity_2d(5, 5, vel);
    require(near(vort[2 * 5 + 2], 2.0, 1e-12), "LBM vorticity field");
    auto slice = lbm::horizontal_slice(5, 5, rho, 2);
    require(slice.size() == 5 && near(slice[3], 5.0), "LBM slice extraction");
    auto line = lbm::streamline_2d(5, 5, vel, {2, 2}, 0.1, 10);
    require(line.size() > 2, "LBM streamline extraction");
    namespace fs = std::filesystem;
    auto tmp = fs::temp_directory_path() / "cfd_lbm_vis";
    fs::create_directories(tmp);
    lbm::write_vtk_structured_2d(tmp / "field.vtk", 5, 5, rho, vel);
    lbm::write_pgm_scalar(tmp / "field.pgm", 5, 5, rho);
    require(fs::file_size(tmp / "field.vtk") > 100 && fs::file_size(tmp / "field.pgm") > 20,
            "LBM headless VTK/PGM output");
    fs::remove_all(tmp);
}

void test_hpc_extensions() {
    using namespace cfd::core;
    TaskGraph g;
    int value = 0;
    auto a = g.add("a", [&] { value = 2; });
    auto b = g.add("b", [&] { value *= 3; }, {a});
    g.add("c", [&] { value += 4; }, {b});
    g.run(true);
    require(value == 10, "task graph dependency execution");
    SimdPack<double, 4> p(2.0), q(3.0);
    require(near((p * q + p).sum(), 32.0), "SIMD pack algebra");
    ScratchArena arena(256);
    auto x = arena.allocate<double>(8);
    x[3] = 4.5;
    require(near(x[3], 4.5) && arena.used() >= 64, "scratch arena");
    arena.reset();
    require(arena.used() == 0, "scratch arena reset");
    CsrBuilder cb(4, 4);
    for (std::size_t i = 0; i < 4; ++i) {
        cb.add(i, i, 4);
        if (i)
            cb.add(i, i - 1, -1);
        if (i + 1 < 4)
            cb.add(i, i + 1, -1);
    }
    auto A = cb.build();
    std::vector<double> exact{1, 2, 3, 4}, rhs(4), sol(4);
    A.multiply(exact, rhs);
    auto ref = mixed_precision_iterative_refinement(A, rhs, sol, 10, 80, 1e-10);
    require(ref.converged, "mixed precision refinement convergence");
    for (std::size_t i = 0; i < 4; ++i)
        require(std::abs(sol[i] - exact[i]) < 1e-8, "mixed precision refinement solution");
    const std::size_t n = 31;
    PoissonMultigrid2D mg(n);
    std::vector<double> f(n * n), u(n * n);
    const double pi = 3.14159265358979323846;
    for (std::size_t j = 0; j < n; ++j)
        for (std::size_t i = 0; i < n; ++i) {
            const double xx = double(i + 1) / double(n + 1), yy = double(j + 1) / double(n + 1);
            f[j * n + i] = 2 * pi * pi * std::sin(pi * xx) * std::sin(pi * yy);
        }
    auto mr = mg.solve(f, u, 60, 1e-8);
    require(mr.converged, "geometric multigrid convergence");
    double err = 0, norm = 0;
    for (std::size_t j = 0; j < n; ++j)
        for (std::size_t i = 0; i < n; ++i) {
            double xx = double(i + 1) / double(n + 1), yy = double(j + 1) / double(n + 1),
                   e = std::sin(pi * xx) * std::sin(pi * yy);
            err += (u[j * n + i] - e) * (u[j * n + i] - e);
            norm += e * e;
        }
    require(std::sqrt(err / norm) < 0.01, "geometric multigrid accuracy");
    std::vector<double> touched(4096);
    numa_first_touch<double>(touched, 2.0);
    require(touched.front() == 2.0 && touched.back() == 2.0, "NUMA first-touch baseline");
    AmgPreconditionerAdapter amg([](std::span<const double> r, std::span<double> z) {
        for (std::size_t i = 0; i < r.size(); ++i)
            z[i] = 0.25 * r[i];
    });
    std::vector<double> rr{4, 8}, zz(2);
    amg(rr, zz);
    require(near(zz[0], 1) && near(zz[1], 2), "AMG adapter callback");
    DistributedCsrPartition left(A, 0, 2), right(A, 2, 4);
    std::vector<double> yl(2), yr(2), hl, hr;
    for (auto c : left.halo_columns())
        hl.push_back(exact[c]);
    for (auto c : right.halo_columns())
        hr.push_back(exact[c]);
    left.multiply(std::span<const double>(exact.data(), 2), hl, yl);
    right.multiply(std::span<const double>(exact.data() + 2, 2), hr, yr);
    require(near(yl[0], rhs[0]) && near(yl[1], rhs[1]) && near(yr[0], rhs[2]) && near(yr[1], rhs[3]),
            "distributed halo-aware CSR SpMV");
    auto rm = cfd::fvm::make_cartesian_hexa_mesh(4, 4, 2);
    auto order = cfd::fvm::morton_cell_order(rm);
    require(order.size() == rm.cell_count(), "Morton cell ordering");
}

void test_fdtd_chemistry_extensions() {
    using namespace cfd;
    fdtd::NearFieldSurfaceSample sample;
    sample.normal = {0, 0, 1};
    sample.area = 1.0;
    sample.electric.x = {1.0, 0.0};
    sample.magnetic.y = {1.0 / 376.730313668, 0.0};
    const std::vector<fdtd::NearFieldSurfaceSample> one{sample}, two{sample, sample};
    const auto f1 = fdtd::nf2ff(one, {0, 0, 1}, 1.0e9);
    const auto f2 = fdtd::nf2ff(two, {0, 0, 1}, 1.0e9);
    require(f1.radiation_intensity > 0.0 && near(f2.radiation_intensity, 4.0 * f1.radiation_intensity, 1e-12),
            "NF2FF area scaling");
    fdtd::SarGrid3D sar;
    sar.nx = sar.ny = sar.nz = 3;
    sar.dx = sar.dy = sar.dz = 0.01;
    const std::size_t cells = 27;
    sar.ex.assign(cells, 10.0);
    sar.ey.assign(cells, 0.0);
    sar.ez.assign(cells, 0.0);
    sar.conductivity.assign(cells, 1.0);
    sar.density.assign(cells, 1000.0);
    const auto local = fdtd::local_sar(sar), s1 = fdtd::sar_1g(sar), s10 = fdtd::sar_10g(sar);
    for (std::size_t i = 0; i < cells; ++i) {
        require(near(local[i], 0.1, 1e-12), "local SAR");
        require(near(s1[i], 0.1, 1e-12), "1g SAR");
        require(near(s10[i], 0.1, 1e-12), "10g SAR");
    }

    using namespace chemistry;
    ReactionNetwork ordered({{"A", 1, 0}, {"B", 1, 0}, {"M", 1, 0}});
    ElementaryReaction ro{{{0, 1}}, {{1, 1}}, {2.0, 0.0, 0.0}, {}};
    ro.forward_orders = {{0, 0.5}};
    ordered.add_reaction(ro);
    auto rates = ordered.reaction_rates(std::vector<double>{4.0, 0.0, 0.0}, 300.0);
    require(near(rates[0], 4.0, 1e-12), "general reaction order");
    ReactionNetwork tb({{"A", 1, 0}, {"B", 1, 0}, {"M", 1, 0}});
    ElementaryReaction tr{{{0, 1}}, {{1, 1}}, {2.0, 0.0, 0.0}, {}};
    tr.third_body = true;
    tr.third_body_efficiencies = {{2, 2.0}};
    tb.add_reaction(tr);
    const double rlow = tb.reaction_rates(std::vector<double>{1.0, 0.0, 1.0}, 300.0)[0];
    const double rhigh = tb.reaction_rates(std::vector<double>{1.0, 0.0, 3.0}, 300.0)[0];
    require(near(rlow, 6.0, 1e-12) && near(rhigh, 14.0, 1e-12), "third-body efficiency");
    ReactionNetwork fall({{"A", 1, 0}, {"B", 1, 0}, {"M", 1, 0}});
    ElementaryReaction fr{{{0, 1}}, {{1, 1}}, {100.0, 0.0, 0.0}, {}};
    fr.low_pressure_limit = ArrheniusRate{10.0, 0.0, 0.0};
    fall.add_reaction(fr);
    const std::vector<double> cf{1.0, 0.0, 9.0};
    const double expected = 100.0 * (10.0 * 10.0 / 100.0) / (1.0 + 10.0 * 10.0 / 100.0);
    require(near(fall.reaction_rates(cf, 300.0)[0], expected, 1e-12), "Lindemann falloff");

    IdealSolutionPhase sol({{"A", 0.0}, {"B", 0.0}});
    const auto mu = sol.chemical_potentials(std::vector<double>{0.25, 0.75}, 300.0, 101325.0);
    require(near(mu[1] - mu[0], gas_constant * 300.0 * std::log(3.0), 1e-12), "ideal solution chemical potential");
    IdealGasPhase gas({{"A", 0.0}, {"B", 0.0}});
    const auto mug = gas.chemical_potentials(std::vector<double>{0.5, 0.5}, 300.0, 202650.0);
    require(near(mug[0], gas_constant * 300.0 * std::log(1.0), 1e-12), "ideal gas activity");

    ReactionNetwork inert({{"A", 1.0, 0}});
    WellStirredReactorNetwork net(inert);
    net.add_reactor({1.0, 300.0, {1.0}});
    net.add_reactor({1.0, 300.0, {0.0}});
    net.connect({0, 1, 0.5});
    const auto amount0 = net.total_species_amounts();
    net.advance(0.5, 0.01);
    const auto amount1 = net.total_species_amounts();
    require(near(amount0[0], amount1[0], 1e-12), "reactor network conservation");
    require(net.reactors()[0].concentrations[0] < 1.0 && net.reactors()[1].concentrations[0] > 0.0,
            "reactor network mixing");

    Nasa7Polynomial nasa{{3.5, 0, 0, 0, 0, 0, 1.0}, 200.0, 5000.0};
    require(near(nasa.cp_molar(500.0), 3.5 * gas_constant, 1e-12), "NASA7 heat capacity");
    require(std::isfinite(nasa.h_molar(500.0)) && std::isfinite(nasa.s_molar(500.0)), "NASA7 thermodynamics");
    DaviesActivity davies;
    const std::vector<int> charge{1, -1};
    const std::vector<double> molality{0.1, 0.1};
    auto gamma = davies.activity_coefficients(charge, molality, 298.15);
    require(gamma[0] < 1.0 && near(gamma[0], gamma[1], 1e-12), "Davies activity");
    SitActivity sit({{0, 1, 0.1}});
    auto sg = sit.activity_coefficients(charge, molality, 298.15);
    require(sg[0] > 0 && sg[1] > 0 && !near(sg[0], gamma[0], 1e-6), "SIT activity correction");
    auto gibbs = minimize_binary_ideal_gibbs(0.0, gas_constant * 300.0 * std::log(3.0), 300.0);
    require(near(gibbs.mole_fraction_a, 0.75, 1e-12) && near(gibbs.mole_fraction_b, 0.25, 1e-12),
            "binary ideal Gibbs minimization");
    require(near(langmuir_surface_coverage(2.0, 0.5), 0.5, 1e-12), "surface complexation Langmuir baseline");
    require(near(binary_ion_exchange_fraction(2.0, 1.0, 0.5), 0.5, 1e-12), "binary ion exchange");

    using namespace electrochemistry;
    PassivationModel pm{0.0, 1.0, 0.2, 0.8, 0.01, 2.0, 1.0};
    require(near(passivation_current_density(0.5, pm), 0.01, 1e-12), "passive plateau");
    require(near(passivation_current_density(0.9, pm), 0.21, 1e-12), "transpassive branch");
    const double recession = faradaic_recession_distance(10.0, 100.0, 0.055845, 7870.0, 2.0);
    require(recession > 0.0, "Faradaic recession");
    ProductLayerState layer{0.0, 100.0, 0.5};
    layer.advance_from_anodic_current(10.0, 100.0, 0.055845, 7870.0, 2.0);
    require(near(layer.thickness, 0.5 * recession, 1e-12) && layer.area_specific_resistance() > 0.0,
            "product-layer growth");
    require(cathodic_protection_current_density(-0.8, -0.5, 0.1) < 0.0, "cathodic protection current");
}

void test_analysis_chemistry_workflow() {
    using namespace cfd;
    std::vector<optics::WavefrontSample> wave;
    for (int ir = 1; ir <= 5; ++ir)
        for (int it = 0; it < 12; ++it) {
            double rho = 0.18 * ir, theta = 2.0 * 3.14159265358979323846 * it / 12.0;
            wave.push_back({rho, theta, 2.0e-6 + 3.0e-6 * optics::zernike(2, 0, rho, theta)});
        }
    auto z = optics::fit_zernike(wave, 2);
    double piston = 0, defocus = 0;
    for (auto c : z) {
        if (c.n == 0 && c.m == 0)
            piston = c.value;
        if (c.n == 2 && c.m == 0)
            defocus = c.value;
    }
    require(std::abs(piston - 2e-6) < 1e-12 && std::abs(defocus - 3e-6) < 1e-12, "Zernike fit");
    require(optics::wavefront_rms(wave) > 0 &&
                near(optics::longitudinal_chromatic_shift(std::vector<double>{0.1, 0.101, 0.099}), 0.002, 1e-12),
            "wavefront/chromatic analysis");
    const double no = 1.5, ne = 1.6;
    require(near(optics::uniaxial_extraordinary_index(no, ne, 0.0), no, 1e-12), "uniaxial optic-axis limit");
    const double retard = optics::birefringent_retardance(no, ne, 1e-3, 500e-9);
    auto j = optics::apply_linear_retarder({{0, 0}, {1, 0}}, retard, 0);
    require(near(std::abs(j.y), 1.0, 1e-12), "birefringent retardance power");
    optics::NonSequentialScene scene;
    scene.add_plane({{0, 0, 0}, {0, 0, -1}, 1, 1, 0.2, 0.5, 0.2, 0.1});
    auto paths = scene.trace({{{0, 0, -1}, {0, 0, 1}, 550}, 1, 0}, 1, 1e-12);
    double power = 0;
    for (auto p : paths)
        power += p.power;
    require(near(power, 0.8, 1e-12), "nonsequential absorption accounting");
    require(paths.size() == 3, "nonsequential scatter branch");
    auto stray = optics::summarize_stray_light(paths, 1.0);
    require(near(stray.total_power, 0.8, 1e-12) && near(stray.absorbed_power, 0.2, 1e-12),
            "stray-light power accounting");

    auto eq = chemistry::equilibrate_binary_salt(2.0, 2.0, 0.0, 1.0);
    require(near(eq.dissolved_a, 1.0, 1e-12) && near(eq.solid_amount, 1.0, 1e-12) &&
                near(eq.saturation_ratio, 1.0, 1e-12),
            "binary precipitation equilibrium");
    auto diss = chemistry::equilibrate_binary_salt(0.1, 0.1, 1.0, 1.0);
    require(diss.reaction_extent < 0 && diss.solid_amount < 1.0, "binary dissolution equilibrium");
    chemistry::ReactionNetwork rp({{"A", 1, 0}, {"B", 1, 0}});
    rp.add_reaction({{{0, 1}}, {{1, 1}}, {2.0, 0, 0}, {}});
    auto edges = chemistry::reaction_path(rp, std::vector<double>{3.0, 0.0}, 300.0);
    require(edges.size() == 1 && edges[0].from == 0 && edges[0].to == 1 && near(edges[0].molar_flux, 6.0, 1e-12),
            "reaction-path flux");
    namespace fs = std::filesystem;
    auto csv = fs::temp_directory_path() / "cfd_thermo.csv";
    {
        std::ofstream o(csv);
        o << "A,-1000\nB,-2000\n";
    }
    auto db = chemistry::load_thermo_csv(csv);
    fs::remove(csv);
    require(db.size() == 2 && near(db[1].standard_chemical_potential, -2000, 1e-12), "thermochemical CSV import");

    auto sweep =
        workflow::parameter_sweep(std::vector<double>{-1, 0, 1, 2}, [](double x) { return (x - 1) * (x - 1); });
    require(sweep[2].objective == 0, "parameter sweep");
    auto inv = workflow::inverse_coordinate_search(
        [](std::span<const double> x) { return (x[0] - 3) * (x[0] - 3) + (x[1] + 2) * (x[1] + 2); }, {0, 0}, {1, 1},
        200, 1e-8);
    require(inv.converged && std::abs(inv.parameters[0] - 3) < 1e-6 && std::abs(inv.parameters[1] + 2) < 1e-6,
            "inverse runner");
    auto checkpoint = fs::temp_directory_path() / "cfd_workflow.chk";
    fs::remove(checkpoint);
    int state = 0;
    workflow::RestartableWorkflow flow;
    auto a = flow.add("a", [&] { state = 2; });
    flow.add("b", [&] { state *= 5; }, {a});
    flow.run(checkpoint);
    require(state == 10, "workflow first run");
    state = 99;
    flow.run(checkpoint);
    require(state == 99, "workflow restart skips completed nodes");
    fs::remove(checkpoint);
}

void test_cfd_framework_extensions() {
    using namespace cfd::fvm;
    FvScalarMatrix fv(2);
    fv.add_diagonal(0, 2);
    fv.add_coefficient(0, 1, -1);
    fv.add_coefficient(1, 0, -1);
    fv.add_diagonal(1, 2);
    fv.add_source(0, 1);
    auto A = fv.matrix();
    std::vector<double> x{1, 2}, y(2);
    A.multiply(x, y);
    require(near(y[0], 0) && near(y[1], 3) && near(fv.source()[0], 1), "generic FVM matrix");
    auto mesh = make_cartesian_hexa_mesh(12, 10, 1);
    std::vector<Vec3> shear(mesh.cell_count());
    for (std::size_t c = 0; c < shear.size(); ++c) {
        auto p = mesh.cells()[c].center;
        shear[c] = {p.y, 0, 0};
    }
    auto lam = eddy_viscosity(mesh, shear, {TurbulenceModel::laminar});
    auto mix = eddy_viscosity(mesh, shear, {TurbulenceModel::mixing_length, 0.05});
    auto smag = eddy_viscosity(mesh, shear, {TurbulenceModel::smagorinsky, 0.01, 0.17});
    double ml = 0, ms = 0;
    for (double v : mix)
        ml += v;
    for (double v : smag)
        ms += v;
    require(lam[0] == 0 && ml > 0 && ms > 0, "turbulence closures");
    require(parse_turbulence_model("WALE") == TurbulenceModel::wale, "turbulence runtime model");
    VofTransport vof(mesh, {1e-3, 1.0});
    vof.initialize([](Vec3 p) { return p.x < 0.5 ? 1.0 : 0.0; });
    const double v0 = vof.volume();
    vof.set_face_flux(std::vector<double>(mesh.face_count(), 0.0));
    vof.step();
    require(near(vof.volume(), v0, 1e-12), "VOF closed conservation");
    for (double a : vof.alpha())
        require(a >= 0 && a <= 1, "VOF boundedness");
    std::vector<double> planar(mesh.cell_count());
    for (std::size_t c = 0; c < planar.size(); ++c)
        planar[c] = mesh.cells()[c].center.x;
    auto csf = csf_surface_tension(mesh, planar, 0.072);
    double middle = 0;
    std::size_t count = 0;
    for (std::size_t c = 0; c < csf.size(); ++c) {
        auto p = mesh.cells()[c].center;
        if (p.x > 0.2 && p.x < 0.8 && p.y > 0.2 && p.y < 0.8) {
            middle += magnitude(csf[c]);
            ++count;
        }
    }
    require(middle / static_cast<double>(count) < 1e-8, "planar CSF curvature");
    auto n = contact_angle_normal({1, 0, 0}, {0, 1, 0}, 3.14159265358979323846 / 2);
    require(std::abs(dot(n, {0, 1, 0})) < 1e-12, "contact-angle normal");
    ParticleCloud cloud;
    cloud.add({{0.5, 0.5, 0.5}, {0, 0, 0}, 1e-3, 1000});
    std::vector<Vec3> fluid(mesh.cell_count(), {1, 0, 0});
    cloud.advance(mesh, fluid, 1e-3, 0.01, true);
    require(cloud.particles()[0].velocity.x > 0 && cloud.particles()[0].position.x > 0.5, "particle drag/tracking");
    double reaction = 0;
    for (auto q : cloud.fluid_momentum_source())
        reaction += q.x * mesh.cells()[0].volume;
    require(reaction < 0, "two-way particle coupling");
}

void test_fem_fdtd_material_extensions() {
    using namespace cfd;
    fem::LinearElastic1D elastic(200.0);
    require(near(elastic.update(0.01).stress, 2.0), "material-law interface");
    fem::BilinearPlasticity1D plastic(200.0, 1.0, 20.0);
    auto p = plastic.update(0.02);
    require(p.equivalent_plastic_strain > 0 && p.tangent < 200.0, "plasticity return mapping");
    auto nh = fem::compressible_neo_hookean({1, 1, 1}, 10, 20);
    require(near(nh.energy_density, 0, 1e-12) && near(nh.nominal_stress[0], 0, 1e-12), "Neo-Hookean reference state");
    auto nh2 = fem::compressible_neo_hookean({1.1, 1, 1}, 10, 20);
    require(nh2.energy_density > 0 && nh2.nominal_stress[0] > 0, "hyperelastic response");
    require(near(fem::penalty_contact_pressure(-1e-3, 1e6), 1000.0), "penalty contact");
    fem::PhaseChangeMaterial pcm{300, 310, 1000, 1200, 200000};
    double h = pcm.specific_enthalpy(305);
    require(pcm.liquid_fraction(305) > 0.49 && pcm.liquid_fraction(305) < 0.51, "phase fraction");
    require(std::abs(pcm.temperature_from_enthalpy(h) - 305) < 1e-8, "enthalpy inversion");
    require(pcm.effective_heat_capacity(305) > 10000, "latent heat capacity");
    require(near(fem::porous_effective_conductivity(.25, 4, 1), 3.25), "porous effective heat transport");
    auto mesh = fem::make_rectangle_tri_mesh(14, 14);
    fem::EddyCurrent2D eddy(mesh, {60, 5000, 60, 1e-10});
    const double pi = 3.14159265358979323846, nu = 1.0, sigma = .1, w = 2 * pi * 60;
    eddy.solve([&](fem::Node2) { return nu; }, [&](fem::Node2) { return sigma; },
               [&](fem::Node2 x) {
                   double a = std::sin(pi * x.x) * std::sin(pi * x.y);
                   return std::complex<double>((2 * pi * pi * nu) * a, w * sigma * a);
               });
    double err = 0, norm = 0;
    for (std::size_t i = 0; i < mesh.node_count(); ++i) {
        double exact = std::sin(pi * mesh.nodes[i].x) * std::sin(pi * mesh.nodes[i].y);
        err += std::norm(eddy.vector_potential()[i] - std::complex<double>(exact, 0));
        norm += exact * exact;
    }
    require(std::sqrt(err / norm) < 0.03, "eddy-current harmonic FEM");
    physics::Material mat{"anisotropic dielectric", 1000, 1000, 1, 0, {2, 3, 4}, {1, 1, 1}, 1.5};
    physics::validate_material(mat);
    fdtd::Maxwell3DConfig mc;
    mc.nx = mc.ny = mc.nz = 8;
    mc.epsilon_rx = 2;
    mc.epsilon_ry = 3;
    mc.epsilon_rz = 4;
    fdtd::Maxwell3D fd(mc);
    fd.initialize_gaussian_ez();
    double e0 = fd.energy();
    fd.step(3);
    require(fd.energy() > 0 && std::isfinite(fd.energy()) && e0 > 0, "anisotropic FDTD material");
}

void test_coupled_cfd_extensions() {
    using namespace cfd;
    chemistry::ReactionNetwork net({{"A", 1, 0}, {"B", 1, 0}});
    net.add_reaction({{{0, 1}}, {{1, 1}}, {2.0, 0, 0}, {}});
    std::vector<std::vector<double>> states{{3, 0}, {1, 0}};
    std::vector<double> temp{300, 300};
    auto src = fvm::reacting_cell_sources(net, states, temp);
    require(near(src[0][0], -6) && near(src[0][1], 6) && near(src[1][0], -2), "reacting CFD source coupling");
    auto cht = multiphysics::conjugate_interface_heat_flux(400, 10, .01, 300, 5, .02);
    require(near(cht.heat_flux, 20000.0, 1e-12) && near(cht.interface_temperature, 380.0, 1e-12),
            "conjugate heat interface");
    require(near(fvm::gray_surface_exchange(300, .8, 300, .7), 0.0, 1e-12), "gray radiation equilibrium");
    require(fvm::gray_surface_exchange(500, .8, 300, .8) > 0, "gray radiation direction");
    require(fvm::optically_thin_radiation_source(500, 300, 1.0) < 0, "participating radiation cooling");
    auto mesh = fvm::make_cartesian_hexa_mesh(3, 2, 1);
    auto moved = fvm::rigidly_moved_mesh(mesh, {1, 2, 3}, {0, 0, 1}, 0.2);
    auto phi = fvm::ale_mesh_flux(mesh, moved, 0.5);
    require(phi.size() == mesh.face_count(), "ALE mesh flux");
    std::vector<double> field(mesh.cell_count());
    for (std::size_t c = 0; c < field.size(); ++c)
        field[c] = 1.0 + mesh.cells()[c].center.x;
    auto remap = fvm::conservative_nearest_remap(mesh, field, moved);
    double i0 = 0, i1 = 0;
    for (std::size_t c = 0; c < field.size(); ++c)
        i0 += field[c] * mesh.cells()[c].volume;
    for (std::size_t c = 0; c < remap.size(); ++c)
        i1 += remap[c] * moved.cells()[c].volume;
    require(near(i0, i1, 1e-12), "conservative mesh remap");
    auto heat = multiphysics::electrochemical_heat_source(std::vector<double>{10, 20}, std::vector<double>{.1, .2},
                                                          std::vector<double>{1, 2});
    require(near(heat[0], 2) && near(heat[1], 6), "electrochemistry heat coupling");
}

void test_lbm_advanced_extensions() {
    using namespace cfd;
    lbm::RegularizedD2Q9Solver reg({40, 40, .75F});
    reg.initialize_taylor_green(.02F);
    double m0 = reg.mass(), e0 = reg.kinetic_energy();
    reg.step(40);
    require(std::abs(reg.mass() - m0) / m0 < 2e-6 && reg.kinetic_energy() < e0, "regularized LBM collision");
    lbm::MrtD2Q9Solver mrt({40, 40, .75F});
    mrt.initialize_taylor_green(.02F);
    m0 = mrt.mass();
    e0 = mrt.kinetic_energy();
    mrt.step(40);
    require(std::abs(mrt.mass() - m0) / m0 < 2e-6 && mrt.kinetic_energy() < e0, "MRT LBM collision");
    io::TriangleSurface cube;
    using P = fem::Point3;
    P p000{.25, .25, .25}, p100{.75, .25, .25}, p010{.25, .75, .25}, p110{.75, .75, .25}, p001{.25, .25, .75},
        p101{.75, .25, .75}, p011{.25, .75, .75}, p111{.75, .75, .75};
    auto tri = [&](P a, P b, P c) { cube.triangles.push_back({a, b, c}); };
    tri(p000, p110, p100);
    tri(p000, p010, p110);
    tri(p001, p101, p111);
    tri(p001, p111, p011);
    tri(p000, p100, p101);
    tri(p000, p101, p001);
    tri(p010, p011, p111);
    tri(p010, p111, p110);
    tri(p000, p001, p011);
    tri(p000, p011, p010);
    tri(p100, p110, p111);
    tri(p100, p111, p101);
    auto mask = lbm::voxelize_surface(cube, {8, 8, 8, {0, 0, 0}, {1, 1, 1}});
    auto idx = [](std::size_t i, std::size_t j, std::size_t k) { return (k * 8 + j) * 8 + i; };
    require(mask[idx(4, 4, 4)] == 1 && mask[idx(0, 0, 0)] == 0, "triangle voxelization");
    auto moved = lbm::transform_surface(cube, {.1, 0, 0}, {0, 0, 1}, 0.0);
    auto mask2 = lbm::voxelize_surface(moved, {8, 8, 8, {0, 0, 0}, {1, 1, 1}});
    require(mask2 != mask, "moving geometry re-voxelization");
}

void test_cpu_solver_gap_closures() {
    using namespace cfd;
    {
        fvm::Compressible1DConfig cfg;
        cfg.cells = 160;
        cfg.length = 1.0;
        cfg.cfl = 0.4;
        cfg.periodic = true;
        fvm::CompressibleEuler1D euler(cfg);
        euler.initialize({1.0, 0.0, 1.0}, {0.125, 0.0, 0.1}, 0.5);
        const double m0 = euler.mass(), E0 = euler.total_energy();
        euler.run(0.03);
        require(std::abs(euler.mass() - m0) / m0 < 1e-11, "compressible mass conservation");
        require(std::abs(euler.total_energy() - E0) / E0 < 1e-11, "compressible energy conservation");
        double pmin = 1e300, rmin = 1e300;
        for (std::size_t i = 0; i < euler.conserved().size(); ++i) {
            auto q = euler.primitive(i);
            pmin = std::min(pmin, q.pressure);
            rmin = std::min(rmin, q.density);
        }
        require(pmin > 0 && rmin > 0, "compressible shock positivity");
        fvm::IdealGasEquationOfState eos;
        require(near(eos.temperature(1.0, 287.05), 1.0, 1e-12), "ideal-gas EOS");
    }
    require(fvm::spalart_allmaras_eddy_viscosity(1e-5, 5e-5) > 0, "SA viscosity helper");
    require(near(fvm::k_epsilon_eddy_viscosity(1.0, 2.0, 4.0), 0.09, 1e-12), "k-epsilon viscosity helper");
    require(fvm::k_omega_sst_eddy_viscosity(1.0, 1.0, 2.0, 3.0) > 0, "SST viscosity helper");
    require(near(fvm::des_length_scale(1.0, 0.1), 0.065, 1e-12), "DES scale helper");
    auto rs = fvm::boussinesq_reynolds_stress(1.0, 3.0, 0.2, 1.0, -1.0, 0.5);
    require(rs.xx < 2.0 && rs.yy > 2.0 && rs.xy < 0, "Reynolds-stress closure helper");
    auto pc = fvm::linear_mushy_phase_change(305, 300, 310, 1000, 1200, 2e5);
    require(pc.liquid_fraction > 0.49 && pc.liquid_fraction < 0.51 && pc.effective_cp > 10000, "FVM phase change");
    require(fvm::drift_flux_slip_velocity(1000, 900, 1e-3, 1e-3) > 0, "drift-flux slip");
    require(fvm::euler_euler_drag_source(0, 1, 2, 3) > 0, "Euler-Euler drag");
    require(fvm::lubrication_thin_film_flux(1e-3, -100, 1e-3) > 0, "thin-film lubrication");

    {
        auto mesh = fem::make_rectangle_tri_mesh(18, 14);
        fem::PenaltyStokes2D stokes(mesh, {1.0, 1e3, 6000, 1e-10});
        const double pi = 3.14159265358979323846;
        stokes.set_dirichlet([&](fem::Node2 p) { return fem::Displacement2{std::sin(pi * p.y), 0.0}; });
        stokes.solve([&](fem::Node2 p) { return fem::Displacement2{pi * pi * std::sin(pi * p.y), 0.0}; });
        double err = 0, norm = 0;
        for (std::size_t i = 0; i < mesh.node_count(); ++i) {
            double ex = std::sin(pi * mesh.nodes[i].y);
            double d = stokes.velocity()[i].x - ex;
            err += d * d;
            norm += ex * ex;
        }
        require(std::sqrt(err / norm) < 0.04, "FEM penalty Stokes manufactured solution");
        require(stokes.divergence_l2() < 1e-2, "FEM penalty Stokes incompressibility");
        auto tr = fem::green_lagrange_truss(1.0, 1.1, 0.01, 200e9);
        require(tr.green_lagrange_strain > 0 && tr.axial_force > 0 && tr.tangent > 0, "nonlinear geometry truss");
    }
    {
        physics::Material m{"shared", 1000, 1000, 1, 5.0, {2, 3, 4}, {1, 1, 1}, 1.6};
        fdtd::Maxwell3DConfig c;
        c.nx = c.ny = c.nz = 8;
        auto mapped = physics::apply_to_fdtd(c, m);
        require(mapped.epsilon_rx == 2 && mapped.epsilon_ry == 3 && mapped.epsilon_rz == 4, "material bridge FDTD");
        require(near(physics::fem_conductivity(m), 5.0) && near(physics::optical_refractive_index(m), 1.6),
                "shared material bridge FEM optics");
    }
    {
        fvm::Particle a{{0, 0, 0}, {1, 0, 0}, 1e-3, 1000}, b{{1e-3, 0, 0}, {-1, 0, 0}, 1e-3, 1000};
        fvm::collide_particles_elastic(a, b, 1.0);
        require(a.velocity.x < 0 && b.velocity.x > 0, "particle collision");
        auto merged = fvm::coalesce_particles(a, b);
        require(merged.diameter > 1e-3, "particle coalescence");
        require(fvm::breakup_by_weber(1000, 10, 1e-3, 0.072), "particle breakup criterion");
        require(fvm::evaporate_d2_law(1e-3, 1e-8, 0.01) < 1e-3, "droplet evaporation");
        auto lift = fvm::saffman_lift_force({1, 0, 0}, {0, 0, 10}, 1000, 1e-6, 1e-3);
        require(lift.y < 0, "particle lift force");
        auto vm = fvm::virtual_mass_force({1, 0, 0}, {0, 0, 0}, 1000, 1e-3);
        require(vm.x > 0, "particle virtual mass");
        auto mesh = fvm::make_cartesian_hexa_mesh(2, 2, 1);
        auto alpha = fvm::particle_volume_fraction(mesh, std::vector<fvm::Particle>{merged});
        double sum = 0;
        for (double v : alpha)
            sum += v;
        require(sum > 0, "dense particle volume fraction");
    }
}

void test_next_phase_portable() {
    using namespace cfd;
    {
        multiphysics::BlockCoupledLinearSystem system({1, 1});
        system.add(0, 0, 0, 0, 2.0);
        system.add(0, 0, 1, 0, 1.0);
        system.add(1, 0, 0, 0, 1.0);
        system.add(1, 0, 1, 0, 3.0);
        system.add_rhs(0, 0, 1.0);
        system.add_rhs(1, 0, 2.0);
        std::vector<double> x(2, 0.0);
        auto result = system.solve(x);
        require(result.converged && near(x[0], 0.2, 1e-10) && near(x[1], 0.6, 1e-10),
                "monolithic block coupling solve");
    }
    {
        optics::SequentialOpticalSystem system(1.0);
        optics::SequentialSurface a;
        a.vertex_z = 0.2;
        a.aperture_radius = 1.0;
        a.refractive_index_after = 1.5;
        optics::SequentialSurface b;
        b.vertex_z = 0.5;
        b.aperture_radius = 1.0;
        b.refractive_index_after = 1.0;
        system.add_surface(a);
        system.add_surface(b);
        auto deformed = multiphysics::deform_optical_system_axially(system, [](double z) { return 1.0e-3 * z; });
        require(near(deformed.surfaces()[0].vertex_z, 0.2002, 1e-12) &&
                    near(deformed.surfaces()[1].vertex_z, 0.5005, 1e-12),
                "deformation to optical surfaces");
    }
    {
        std::vector<fvm::IdealGasSpeciesThermo> sp{{0.028, 1000.0, 0.0}, {0.032, 900.0, -2.0e5}};
        std::vector<double> y{0.75, 0.25};
        auto st = fvm::ideal_gas_mixture_state(1200.0, 101325.0, y, sp);
        require(st.density > 0 && st.gamma > 1.0 && st.cp > st.cv, "combustion ideal-gas thermo state");
        require(near(fvm::temperature_from_mixture_enthalpy(st.specific_enthalpy, y, sp), 1200.0, 1e-10),
                "combustion enthalpy inversion");
        fvm::ResidualHistory history;
        fvm::FunctionObjectManager manager;
        manager.add([&](const fvm::SolverReport& r) { history.observe(r); });
        manager.execute({5, 0.1, {{"p", 1e-6}, {"U", 2e-6}}});
        require(history.records().size() == 1 && near(history.latest("p"), 1e-6), "residual/function object framework");
    }
    {
        fdtd::CylindricalTM solver({36, 72, 1e-3, 1e-3, 0.25, 1.0, 1.0});
        solver.initialize_gaussian_ez(0.0, 0.0355, 0.006, 1.0);
        const double e0 = solver.energy();
        solver.step(80);
        const double e1 = solver.energy();
        require(e0 > 0 && e1 > 0 && std::isfinite(e1) && std::isfinite(solver.max_field()),
                "cylindrical TM FDTD finite energy");
        require(e1 / e0 < 5.0 && solver.max_field() < 5.0, "cylindrical TM FDTD stability");
        const auto& ez = solver.ez();
        const std::size_t nr = 36, nz = 72;
        double symmetry = 0.0;
        for (std::size_t j = 0; j < nz - 1; ++j) {
            std::size_t mirror = (nz - 2) - j;
            symmetry = std::max(symmetry, std::abs(ez[j * nr] - ez[mirror * nr]));
        }
        require(symmetry < 2e-10, "cylindrical TM axial symmetry regression");
    }
}

void test_portable_physics_depth() {
    using namespace cfd;
    {
        const double dt = 0.1, t = 1.0;
        require(near(fvm::bdf2_derivative(t * t, (t - dt) * (t - dt), (t - 2 * dt) * (t - 2 * dt), dt), 2 * t, 1e-12),
                "BDF2 quadratic derivative");
        require(near(fvm::crank_nicolson_update(1.0, 2.0, 4.0, 0.5), 2.5, 1e-12), "Crank-Nicolson trapezoid update");
    }
    {
        std::vector<double> T{500, 300}, eps{0.8, 0.7}, area{1, 1}, vf{0, 1, 1, 0};
        auto r = fvm::solve_gray_radiosity(T, eps, area, vf);
        double expected = fvm::gray_surface_exchange(500, 0.8, 300, 0.7);
        require(near(r.net_heat[0], expected, 1e-10 * std::abs(expected)) &&
                    near(r.net_heat[0] + r.net_heat[1], 0.0, 1e-8),
                "gray radiosity enclosure conservation");
        fvm::OpticallyThinRadiation model(0.2);
        require(model.volumetric_source(500, 300) < 0, "radiation model interface");
    }
    {
        fvm::PremixedFlame1DConfig cfg;
        cfg.cells = 301;
        cfg.length = 1.0;
        cfg.diffusivity = 1e-3;
        cfg.reaction_rate = 1.0;
        cfg.dt = 1e-3;
        fvm::PremixedFlame1D flame(cfg);
        flame.initialize_front(0.25, 0.025);
        double x0 = flame.front_location();
        flame.step(1500);
        double x1 = flame.front_location();
        require(x1 > x0 + 0.03 && x1 < 0.5, "premixed flame front propagation");
        for (double c : flame.progress())
            require(c >= 0 && c <= 1, "premixed flame bounded progress");
        require(flame.temperature(0) > flame.temperature(flame.progress().size() - 1), "premixed flame thermo mapping");
    }
    {
        chemistry::ShomatePolynomial n2{28.98641,  1.853978, -9.647459, 16.63537, 0.000117,
                                        -8.671914, 226.4168, 0.0,       298.0,    1000.0};
        require(std::abs(n2.cp_molar(298.15) - 29.1248) < 0.01, "Shomate N2 heat capacity");
        require(std::isfinite(n2.h_molar_kj(500)) && std::isfinite(n2.s_molar(500)), "Shomate enthalpy entropy");
        chemistry::BinaryPitzer11 pitzer{0.0765, 0.2664, 0.00127};
        double g1 = pitzer.mean_activity_coefficient(1e-8), g2 = pitzer.mean_activity_coefficient(1.0);
        require(std::abs(g1 - 1.0) < 1e-3 && g2 > 0 && std::isfinite(g2), "binary Pitzer activity baseline");
    }
    {
        double ip = electrochemistry::marcus_hush_chidsey_current_density(2.0, 0.05, 50000.0);
        double im = electrochemistry::marcus_hush_chidsey_current_density(2.0, -0.05, 50000.0);
        require(near(electrochemistry::marcus_hush_chidsey_current_density(2.0, 0.0, 50000.0), 0.0, 1e-14) &&
                    near(ip, -im, 1e-10 * std::max(1.0, std::abs(ip))),
                "Marcus-Hush-Chidsey symmetry");
        electrochemistry::PhaseFieldCorrosion1D pf(101, 0.01, 0.05, 0.002);
        pf.initialize_interface(0.5, 0.04);
        std::vector<double> drive(101, 0.1);
        auto before = pf.order_parameter();
        pf.step(1e-3, drive);
        double change = 0;
        for (std::size_t i = 0; i < before.size(); ++i) {
            change = std::max(change, std::abs(before[i] - pf.order_parameter()[i]));
            require(pf.order_parameter()[i] >= -1 && pf.order_parameter()[i] <= 1, "phase-field boundedness");
        }
        require(change > 0, "spatial phase-field corrosion evolution");
    }
    {
        std::vector<optics::NsRay> rays{{{{0, 0, 0}, {0, 0, 1}, 550}, 0.1, 3},
                                        {{{0, 0, 0}, {0, 0, 1}, 550}, 0.3, 2},
                                        {{{0, 0, 0}, {0, 0, 1}, 550}, 0.5, 1}};
        auto ghosts = optics::rank_ghost_paths(rays);
        require(ghosts.size() == 2 && near(ghosts[0].power, 0.3) && ghosts[0].bounces == 2, "ghost-path ranking");
    }
}

void test_rf_and_spice() {
    using namespace cfd;
    {
        rf::Matrix2C z{{75.0, 5.0}, {10.0, -2.0}, {10.0, -2.0}, {60.0, 1.0}};
        auto s = rf::z_to_s(z, 50.0);
        auto zr = rf::s_to_z(s, 50.0);
        require(std::abs(zr.a11 - z.a11) < 1e-10 && std::abs(zr.a12 - z.a12) < 1e-10 &&
                    std::abs(zr.a21 - z.a21) < 1e-10 && std::abs(zr.a22 - z.a22) < 1e-10,
                "RF Z/S roundtrip");
        rf::TransmissionLine line;
        line.characteristic_impedance = 50.0;
        line.propagation_constant = {0.0, 2.0 * std::numbers::pi};
        line.length = 0.25;
        auto ls = line.s_parameters(50.0);
        require(std::abs(ls.a11) < 1e-12 && std::abs(std::abs(ls.a21) - 1.0) < 1e-12,
                "matched lossless transmission line");
        auto micro = rf::microstrip_quasi_static(2.0e-3, 1.0e-3, 4.4, 2.4e9);
        require(micro.characteristic_impedance > 40 && micro.characteristic_impedance < 60 &&
                    micro.guided_wavelength > 0,
                "microstrip quasi-static model");
        auto coax = rf::coax_quasi_static(1.0e-3, 2.3e-3, 1.0, 1.0e9);
        require(std::abs(coax.characteristic_impedance - 49.94) < 0.2 && coax.propagation_constant_rad_per_m > 0,
                "coax TEM quasi-static model");
        auto wr90 = rf::rectangular_waveguide_te10(22.86e-3, 10.16e-3, 1.0, 1.0, 10.0e9);
        require(std::abs(wr90.cutoff_frequency_hz / 1.0e9 - 6.557) < 0.01 && wr90.wave_impedance > 376.0,
                "rectangular waveguide TE10 model");
        auto hz = rf::hertzian_dipole_far_field(0.01, 1.0, 1.0e9, std::numbers::pi / 2.0, 10.0);
        require(std::abs(std::abs(hz.e_theta / hz.h_phi) - 376.730313668) < 1e-9, "Hertzian dipole E/H ratio");
        auto hzaxis = rf::hertzian_dipole_far_field(0.01, 1.0, 1.0e9, 0.0, 10.0);
        require(std::abs(hzaxis.e_theta) < 1e-15 && near(rf::shielding_effectiveness_db(1.0, 0.1), 20.0, 1e-12),
                "Hertzian axis/shielding utility");
        auto strip = rf::stripline_quasi_static(0.5e-3, 1.0e-3, 4.0, 2.4e9);
        require(strip.characteristic_impedance > 45 && strip.characteristic_impedance < 55 &&
                    near(strip.effective_permittivity, 4.0),
                "stripline quasi-static model");
        auto cpw = rf::coplanar_waveguide_quasi_static(1.0e-3, 0.2e-3, 4.0, 2.4e9);
        require(cpw.characteristic_impedance > 20 && cpw.characteristic_impedance < 100 &&
                    cpw.effective_permittivity > 1.0,
                "CPW quasi-static model");
        std::istringstream touch{"# GHZ S RI R 50\n1 0 0 1 0 1 0 0 0\n"};
        double z0 = 0;
        auto pts = rf::read_touchstone_s2p(touch, &z0);
        require(pts.size() == 1 && near(pts[0].frequency_hz, 1e9) && near(z0, 50) && near(std::abs(pts[0].s.a21), 1),
                "Touchstone S2P reader");
        rf::SinusoidalDipoleSolver dipole(0.5 * 299792458.0 / 1e9, 1e9);
        double rr = dipole.radiation_resistance(), d = dipole.directivity();
        require(std::abs(rr - 73.13) < 0.5, "half-wave dipole radiation resistance");
        require(std::abs(d - 1.64) < 0.03, "half-wave dipole directivity");
        rf::ComplexMatrix zn(3);
        for (std::size_t i = 0; i < 3; ++i)
            zn(i, i) = 50.0;
        zn(0, 1) = zn(1, 0) = 5.0;
        zn(1, 2) = zn(2, 1) = {2.0, 1.0};
        const std::vector<double> zref{50.0, 60.0, 75.0};
        auto sn = rf::z_to_s(zn, zref);
        auto zn2 = rf::s_to_z(sn, zref);
        for (std::size_t r = 0; r < 3; ++r)
            for (std::size_t c = 0; c < 3; ++c)
                require(std::abs(zn2(r, c) - zn(r, c)) < 1e-9, "RF N-port Z/S roundtrip");
        const std::vector<double> zref2{75.0, 50.0, 100.0};
        auto ren = rf::renormalize_s(sn, zref, zref2);
        auto back = rf::renormalize_s(ren, zref2, zref);
        for (std::size_t r = 0; r < 3; ++r)
            for (std::size_t c = 0; c < 3; ++c)
                require(std::abs(back(r, c) - sn(r, c)) < 1e-9, "RF N-port renormalization");
        rf::NPortPoint np;
        np.frequency_hz = 2e9;
        np.s = sn;
        np.reference_impedance = zref;
        std::ostringstream ts;
        rf::write_touchstone(ts, std::span<const rf::NPortPoint>(&np, 1), true);
        std::istringstream tin(ts.str());
        auto nr = rf::read_touchstone(tin);
        require(nr.size() == 1 && nr[0].s.size() == 3 && std::abs(nr[0].s(2, 1) - sn(2, 1)) < 1e-5,
                "Touchstone v2 SnP roundtrip");
        rf::ComplexMatrix line_s(2);
        const rf::Complex gamma{0.1, 12.0};
        const double length = 0.03;
        line_s(0, 1) = line_s(1, 0) = std::exp(-gamma * length);
        const std::array<rf::Complex, 2> gammas{gamma, gamma};
        const std::array<double, 2> shifts{length / 2, length / 2};
        auto deembedded = rf::shift_reference_planes(line_s, gammas, shifts);
        require(std::abs(deembedded(1, 0) - rf::Complex{1, 0}) < 1e-12, "RF reference-plane deembedding");
        require(near(rf::return_loss_db({0.1, 0}), 20.0, 1e-12) && near(rf::vswr({0.1, 0}), 1.1 / 0.9, 1e-12),
                "RF return loss VSWR");
        rf::Matrix2C amp{{0.1, 0}, {0.01, 0}, {2.0, 0}, {0.1, 0}};
        auto stable = rf::stability_metrics(amp);
        require(stable.rollett_k > 1.0 && stable.mu > 1.0, "RF K/mu stability metrics");
        require(near(rf::transducer_gain(amp), 4.0, 1e-12), "RF matched transducer gain");
        auto sc = rf::source_stability_circle(amp), lc = rf::load_stability_circle(amp);
        require(sc.valid && lc.valid && sc.radius >= 0 && lc.radius >= 0, "RF source/load stability circles");
        auto pull = rf::sample_load_pull(amp, {}, 3, 8, 0.8);
        require(pull.size() == 25 && pull.front().transducer_gain > 0, "RF load-pull sampling");
        rf::ComplexMatrix identity4 = rf::ComplexMatrix::identity(4);
        auto mixed = rf::single_ended_to_mixed_mode(identity4);
        for (std::size_t i = 0; i < 4; ++i)
            for (std::size_t j = 0; j < 4; ++j)
                require(std::abs(mixed(i, j) - (i == j ? rf::Complex{1, 0} : rf::Complex{})) < 1e-12,
                        "RF mixed-mode transform");
        auto quality = rf::network_quality(identity4);
        require(quality.passive && quality.reciprocal && std::abs(quality.maximum_singular_value - 1.0) < 1e-12,
                "RF N-port passivity reciprocity");
        rf::ComplexMatrix active2(2);
        active2(1, 0) = 2.0;
        require(!rf::network_quality(active2).passive, "RF active-network passivity rejection");
        const auto match = rf::synthesize_l_match(50, 200);
        require(!match.shunt_first && near(match.series_reactance_ohm, 50 * std::sqrt(3.0), 1e-12) &&
                    near(match.shunt_susceptance_siemens, std::sqrt(3.0) / 200, 1e-12),
                "RF L-match synthesis");
        const double lambda = 299792458.0 / 1.0e9;
        const std::array<rf::ArrayElement, 2> elements{{{0, 0, 0, {1, 0}}, {lambda / 2, 0, 0, {1, 0}}}};
        require(std::abs(rf::array_factor(elements, 1e9, 0, 0) - rf::Complex{1, 0}) < 1e-12 &&
                    std::abs(rf::array_factor(elements, 1e9, std::numbers::pi / 2, 0)) < 1e-12,
                "RF phased-array factor");
        auto pol = rf::polarization_metrics({1, 0}, {0, 1});
        require(std::abs(pol.axial_ratio - 1.0) < 1e-12, "RF circular polarization axial ratio");
        const double mom_lambda = 299792458.0 / 3.0e8;
        rf::ThinWireMomConfig mom_cfg{0.48 * mom_lambda, 0.001 * mom_lambda, 3.0e8, 31U, 1.0, 12U};
        auto mom = rf::solve_center_fed_thin_wire(mom_cfg);
        require(mom.current_a.size() == mom_cfg.segments && std::isfinite(mom.feed_impedance_ohm.real()) &&
                    std::isfinite(mom.feed_impedance_ohm.imag()) && mom.feed_impedance_ohm.real() > 0.0 &&
                    mom.feed_current_a > 0.0,
                "RF thin-wire MoM finite passive feed solution");
        double mom_symmetry = 0.0, mom_peak = 0.0;
        for (std::size_t i = 0; i < mom.current_a.size(); ++i) {
            mom_symmetry =
                std::max(mom_symmetry, std::abs(mom.current_a[i] - mom.current_a[mom.current_a.size() - 1U - i]));
            mom_peak = std::max(mom_peak, std::abs(mom.current_a[i]));
        }
        require(mom_symmetry / std::max(mom_peak, 1e-30) < 1e-10, "RF thin-wire MoM mirror symmetry");
        rf::PeecFilamentSystem peec;
        peec.add_segment({{0, 0, 0}, {0, 0, 0.1}, 5e-4, 5.8e7});
        peec.add_segment({{0.02, 0, 0}, {0.02, 0, 0.1}, 5e-4, 5.8e7});
        auto pm = peec.extract();
        require(pm.size == 2 && pm.resistance_ohm[0] > 0 && pm.partial_inductance_h[0] > 0 &&
                    pm.partial_inductance_h[1] > 0 && pm.partial_inductance_h[1] < pm.partial_inductance_h[0],
                "RF PEEC partial R/L extraction");
        auto pz = peec.impedance_matrix(1e8);
        require(std::abs(pz[1] - pz[2]) < 1e-15, "RF PEEC reciprocal mutual impedance");
        const std::array<rf::Complex, 2> pv{{{1, 0}, {0, 0}}};
        auto pi = peec.solve_currents(0.0, pv);
        require(std::abs(pi[0].real() - 1.0 / pm.resistance_ohm[0]) / (1.0 / pm.resistance_ohm[0]) < 1e-12 &&
                    std::abs(pi[1]) < 1e-12,
                "RF PEEC DC conductor solve");
        const auto& first_segment = peec.segments().front();
        const double skin = rf::skin_depth_m(first_segment.conductivity_s_per_m, 1.0e9);
        const double rac = rf::round_wire_ac_resistance(first_segment, 1.0e9);
        require(skin < first_segment.radius_m && rac > pm.resistance_ohm[0], "RF PEEC skin-effect resistance");
        auto psz = peec.impedance_matrix_skin_effect(1e9);
        require(psz[0].real() > pz[0].real(), "RF PEEC skin-effect impedance matrix");
        rf::PeecFilamentSystem single;
        single.add_segment(first_segment);
        const std::string spice_peec = single.to_spice_subcircuit("WIRE");
        auto peec_circuit = circuit::Circuit::parse_spice(std::string("V1 in 0 1\nXW in 0 WIRE\n") + spice_peec);
        auto peec_op = peec_circuit.dc_operating_point();
        require(peec_op.converged, "RF PEEC SPICE subcircuit export");
    }
    {
        circuit::Circuit c;
        auto vin = c.node("in"), vout = c.node("out");
        c.add_voltage_source("V1", vin, 0, 10.0, {1.0, 0.0});
        c.add_resistor("R1", vin, vout, 1000);
        c.add_resistor("R2", vout, 0, 1000);
        auto op = c.dc_operating_point();
        require(op.converged && near(c.voltage(op, "out"), 5.0, 1e-8), "SPICE MNA DC divider");
        auto ac = c.ac(1000.0, &op);
        require(near(std::abs(c.voltage(ac, "out")), 0.5, 1e-10), "SPICE MNA AC divider");
    }
    {
        circuit::Circuit rc;
        auto in = rc.node("in"), out = rc.node("out");
        rc.add_voltage_source("V", in, 0, 0.0, {1.0, 0.0});
        rc.add_resistor("R", in, out, 1000);
        rc.add_capacitor("C", out, 0, 1e-6);
        const double fc = 1.0 / (2.0 * std::numbers::pi * 1000.0 * 1e-6);
        auto ac = rc.ac(fc);
        auto h = rc.voltage(ac, "out");
        require(std::abs(std::abs(h) - 1 / std::sqrt(2.0)) < 1e-5 &&
                    std::abs(std::arg(h) + std::numbers::pi / 4.0) < 1e-5,
                "SPICE RC low-pass AC");
        auto pole = circuit::estimate_dominant_pole(rc, "out", 1.0, 1.0e5, 101);
        require(pole.found && std::abs(pole.pole_frequency_hz / fc - 1.0) < 0.03, "SPICE dominant pole estimation");
    }
    {
        circuit::Circuit diode;
        auto vin = diode.node("vin"), out = diode.node("out");
        diode.add_voltage_source("V", vin, 0, 5.0);
        diode.add_resistor("R", vin, out, 1000);
        diode.add_diode_model("DM", {1e-12, 1.0, 300});
        diode.add_diode("D", out, 0, "DM");
        auto op = diode.dc_operating_point();
        require(op.converged && diode.voltage(op, "out") > 0.45 && diode.voltage(op, "out") < 0.9,
                "SPICE diode Newton operating point");
        auto temps = diode.dc_temperature_sweep(280.0, 340.0, 3);
        require(temps.size() == 3 && temps.front().operating_point.converged &&
                    temps.back().operating_point.converged &&
                    std::abs(temps.front().operating_point.node_voltage[out] -
                             temps.back().operating_point.node_voltage[out]) > 1e-3,
                "SPICE compact-model temperature sweep");
    }
    {
        const char* net =
            "* divider\nV1 in 0 DC 10 AC 1 0\nR1 in out 1k\nR2 out 0 1k\n.model dm D (IS=1p N=1)\nD1 out 0 dm\n";
        auto c = circuit::Circuit::parse_spice(net);
        auto op = c.dc_operating_point();
        require(op.converged && c.voltage(op, "out") > 0.4 && c.voltage(op, "out") < 1.0, "SPICE netlist/model parser");
        require(near(circuit::parse_spice_number("1meg"), 1e6) && near(circuit::parse_spice_number("10u"), 1e-5),
                "SPICE suffix parser");
    }
    {
        circuit::Circuit mos;
        auto vdd = mos.node("vdd"), gate = mos.node("g"), drain = mos.node("d");
        mos.add_voltage_source("VDD", vdd, 0, 5.0);
        mos.add_voltage_source("VG", gate, 0, 2.0);
        mos.add_resistor("RD", vdd, drain, 2000);
        mos.add_mos_model("N1", {false, 0.7, 1e-3, 0.02});
        mos.add_mosfet("M1", drain, gate, 0, 0, "N1");
        auto op = mos.dc_operating_point();
        require(op.converged && mos.voltage(op, "d") > 0.0 && mos.voltage(op, "d") < 5.0,
                "SPICE level-1 MOS operating point");
    }
    {
        circuit::Circuit rc;
        auto in = rc.node("in"), out = rc.node("out");
        rc.add_voltage_source("V", in, 0, 1.0);
        rc.add_resistor("R", in, out, 1000);
        rc.add_capacitor("C", out, 0, 1e-6);
        auto tr = rc.transient(1e-4, 10);
        require(tr.size() == 11 && tr[1].node_voltage[out] > 0 &&
                    tr.back().node_voltage[out] > tr[1].node_voltage[out] && tr.back().node_voltage[out] < 1.0,
                "SPICE backward-Euler RC transient");
        circuit::Circuit rc2;
        in = rc2.node("in");
        out = rc2.node("out");
        rc2.add_voltage_source("V", in, 0, 1.0);
        rc2.add_resistor("R", in, out, 1000);
        rc2.add_capacitor("C", out, 0, 1e-6);
        const double dt = 2e-4;
        const std::size_t steps = 20;
        const double exact = 1.0 - std::exp(-dt * static_cast<double>(steps) / (1000.0 * 1e-6));
        auto be = rc2.transient(dt, steps, circuit::TransientMethod::backward_euler);
        auto trap = rc2.transient(dt, steps, circuit::TransientMethod::trapezoidal);
        auto bdf = rc2.transient(dt, steps, circuit::TransientMethod::bdf2);
        const double ebe = std::abs(be.back().node_voltage[out] - exact),
                     et = std::abs(trap.back().node_voltage[out] - exact),
                     ebdf = std::abs(bdf.back().node_voltage[out] - exact);
        require(et < ebe && ebdf < ebe, "SPICE trapezoidal/BDF2 transient accuracy");
        circuit::Circuit pulse;
        out = pulse.node("out");
        circuit::PulseWaveform pw{0.0, 2.0, 1e-3, 0.0, 0.0, 1e-3, 3e-3};
        pulse.add_voltage_source("VP", out, 0, 0.0, {}, [pw](double t) { return pw(t); });
        auto ptr = pulse.transient(0.5e-3, 7);
        require(near(ptr[1].node_voltage[out], 0.0, 1e-12) && near(ptr[2].node_voltage[out], 2.0, 1e-12) &&
                    near(ptr[5].node_voltage[out], 0.0, 1e-12),
                "SPICE PULSE transient source");
        const char* wave_net = "VP out 0 PULSE(0 1 1m 0 0 1m 3m)\n";
        auto wp = circuit::Circuit::parse_spice(wave_net);
        auto wpt = wp.transient(0.5e-3, 5);
        const auto wout = wp.find_node("out");
        require(near(wpt[2].node_voltage[wout], 1.0, 1e-12) && near(wpt[5].node_voltage[wout], 0.0, 1e-12),
                "SPICE PULSE netlist parser");
        circuit::SineWaveform sine{1.0, 2.0, 50.0, 0.0, 0.0, 0.0};
        require(near(sine(0.005), 3.0, 1e-12), "SPICE SIN waveform");
        circuit::PwlWaveform pwl{{{0.0, 0.0}, {1.0, 2.0}, {2.0, 0.0}}};
        require(near(pwl(0.5), 1.0, 1e-12) && near(pwl(3.0), 0.0, 1e-12), "SPICE PWL waveform");
    }
    {
        circuit::Circuit dep;
        auto ctrl = dep.node("ctrl"), out = dep.node("out");
        dep.add_voltage_source("VCTRL", ctrl, 0, 1.0);
        dep.add_vcvs("E1", out, 0, ctrl, 0, 2.0);
        dep.add_resistor("RL", out, 0, 1000);
        auto op = dep.dc_operating_point();
        require(op.converged && near(dep.voltage(op, "out"), 2.0, 1e-8), "SPICE VCVS stamp");
        circuit::Circuit gm;
        ctrl = gm.node("ctrl");
        out = gm.node("out");
        gm.add_voltage_source("VC", ctrl, 0, 1.0);
        gm.add_vccs("G1", out, 0, ctrl, 0, 1e-3);
        gm.add_resistor("R", out, 0, 1000);
        op = gm.dc_operating_point();
        require(op.converged && near(gm.voltage(op, "out"), -1.0, 1e-8), "SPICE VCCS stamp");
    }
    {
        circuit::Circuit current;
        auto ctrl = current.node("ctrl"), out = current.node("out");
        current.add_voltage_source("VSENSE", ctrl, 0, 1.0);
        current.add_resistor("RS", ctrl, 0, 1000);
        current.add_cccs("F1", out, 0, "VSENSE", 2.0);
        current.add_resistor("RL", out, 0, 1000);
        auto op = current.dc_operating_point();
        require(op.converged && std::abs(std::abs(current.voltage(op, "out")) - 2.0) < 1e-7, "SPICE CCCS stamp");
        circuit::Circuit trans;
        ctrl = trans.node("ctrl");
        out = trans.node("out");
        trans.add_voltage_source("VSENSE", ctrl, 0, 1.0);
        trans.add_resistor("RS", ctrl, 0, 1000);
        trans.add_ccvs("H1", out, 0, "VSENSE", 1000);
        trans.add_resistor("RL", out, 0, 1000);
        op = trans.dc_operating_point();
        require(op.converged && std::abs(std::abs(trans.voltage(op, "out")) - 1.0) < 1e-7, "SPICE CCVS stamp");
    }
    {
        circuit::Circuit sw;
        auto vin = sw.node("vin"), out = sw.node("out"), ctl = sw.node("ctl");
        sw.add_voltage_source("VIN", vin, 0, 10.0);
        sw.add_voltage_source("VC", ctl, 0, 5.0);
        sw.add_resistor("R", vin, out, 1000);
        sw.add_switch_model("SWM", {1.0, 1e12, 2.5, 0.1});
        sw.add_switch("S1", out, 0, ctl, 0, "SWM");
        auto op = sw.dc_operating_point();
        require(op.converged && sw.voltage(op, "out") < 0.02, "SPICE voltage-controlled switch");
        circuit::Circuit bjt;
        auto vcc = bjt.node("vcc"), base = bjt.node("b"), collector = bjt.node("c");
        bjt.add_voltage_source("VCC", vcc, 0, 5.0);
        bjt.add_voltage_source("VB", base, 0, 0.7);
        bjt.add_resistor("RC", vcc, collector, 1000);
        bjt.add_bjt_model("QN", {false, 1e-15, 100, 1, 300});
        bjt.add_bjt("Q1", collector, base, 0, "QN");
        op = bjt.dc_operating_point_homotopy();
        require(op.converged && bjt.voltage(op, "c") > 0.0 && bjt.voltage(op, "c") < 5.0,
                "SPICE BJT Ebers-Moll baseline");
        circuit::Circuit jf;
        vcc = jf.node("vcc");
        auto drain = jf.node("d");
        jf.add_voltage_source("VDD", vcc, 0, 5.0);
        jf.add_resistor("RD", vcc, drain, 1000);
        jf.add_jfet_model("JN", {false, -2.0, 1e-3, 0.0});
        jf.add_jfet("J1", drain, 0, 0, "JN");
        op = jf.dc_operating_point();
        require(op.converged && jf.voltage(op, "d") > 1.0 && jf.voltage(op, "d") < 5.0, "SPICE JFET baseline");
    }
    {
        circuit::Circuit transformer;
        auto in = transformer.node("in"), out = transformer.node("out");
        transformer.add_voltage_source("V", in, 0, 0.0, {1.0, 0.0});
        transformer.add_inductor("L1", in, 0, 10e-3);
        transformer.add_inductor("L2", out, 0, 10e-3);
        transformer.add_mutual_inductance("K1", "L1", "L2", 0.9);
        transformer.add_resistor("RL", out, 0, 1000);
        auto ac = transformer.ac(1000);
        require(std::abs(transformer.voltage(ac, "out")) > 0.1, "SPICE mutual inductance AC");
        circuit::Circuit rfload;
        auto src = rfload.node("src"), p1 = rfload.node("p1"), p2 = rfload.node("p2");
        rfload.add_voltage_source("VAC", src, 0, 0.0, {1.0, 0.0});
        rfload.add_resistor("RS", src, p1, 50.0);
        rfload.add_resistor("RL", p2, 0, 50.0);
        rf::NPortPoint two;
        two.frequency_hz = 1e9;
        two.s = rf::ComplexMatrix(2);
        two.s(0, 1) = 0.5;
        two.s(1, 0) = 0.5;
        two.reference_impedance = {50, 50};
        rfload.add_nport("SNET", {{p1, 0}, {p2, 0}}, {two});
        auto nac = rfload.ac(1e9);
        require(std::abs(std::abs(rfload.voltage(nac, "p2")) - 0.25) < 1e-10, "RF N-port device in AC MNA");
        circuit::Circuit passive;
        auto a = passive.node("a"), b = passive.node("b");
        passive.add_resistor("RAB", a, b, 100.0);
        passive.add_resistor("RA", a, 0, 200.0);
        passive.add_resistor("RB", b, 0, 200.0);
        auto cs = circuit::two_port_s_parameters(passive, {"a", "0"}, {"b", "0"}, 1e6, 50.0);
        require(std::abs(cs.a12 - cs.a21) < 1e-12 && std::abs(cs.a11) < 1.0 && std::abs(cs.a21) < 1.0,
                "SPICE direct two-port S extraction");
        circuit::Circuit tline;
        auto ta = tline.node("ta"), tb = tline.node("tb");
        const std::array<double, 1> tf{{1.0e9}};
        cfd::rf::TemTransmissionLine tl{50.0, 299792458.0, 0.0, 299792458.0 / (4.0e9)};
        tline.add_tem_transmission_line("TL", {ta, 0}, {tb, 0}, tl, tf, 50.0);
        auto tls = circuit::two_port_s_parameters(tline, {"ta", "0"}, {"tb", "0"}, 1.0e9, 50.0);
        require(std::abs(std::abs(tls.a21) - 1.0) < 1e-10 && std::abs(tls.a11) < 1e-10,
                "SPICE TEM transmission-line N-port element");
        circuit::Circuit mline;
        ta = mline.node("ta");
        tb = mline.node("tb");
        mline.add_microstrip_line("ML", {ta, 0}, {tb, 0}, 2.0e-3, 1.0e-3, 4.4, 0.02, tf, 50.0);
        auto mls = circuit::two_port_s_parameters(mline, {"ta", "0"}, {"tb", "0"}, 1.0e9, 50.0);
        require(std::isfinite(std::abs(mls.a21)) && std::abs(mls.a21) > 0.5, "SPICE microstrip circuit element");
        const auto strip_model = rf::stripline_quasi_static(0.5e-3, 1.0e-3, 4.0, 1.0e9);
        circuit::Circuit sline;
        ta = sline.node("ta");
        tb = sline.node("tb");
        sline.add_stripline_line("SL", {ta, 0}, {tb, 0}, 0.5e-3, 1.0e-3, 4.0, 0.02, tf,
                                 strip_model.characteristic_impedance);
        auto sls = circuit::two_port_s_parameters(sline, {"ta", "0"}, {"tb", "0"}, 1.0e9,
                                                  strip_model.characteristic_impedance);
        require(std::abs(sls.a11) < 1e-9 && std::abs(std::abs(sls.a21) - 1.0) < 1e-9,
                "SPICE stripline circuit element");
        const auto cpw_model = rf::coplanar_waveguide_quasi_static(1.0e-3, 0.2e-3, 4.0, 1.0e9);
        circuit::Circuit cline;
        ta = cline.node("ta");
        tb = cline.node("tb");
        cline.add_coplanar_waveguide_line("CPW", {ta, 0}, {tb, 0}, 1.0e-3, 0.2e-3, 4.0, 0.02, tf,
                                          cpw_model.characteristic_impedance);
        auto cls =
            circuit::two_port_s_parameters(cline, {"ta", "0"}, {"tb", "0"}, 1.0e9, cpw_model.characteristic_impedance);
        require(std::abs(cls.a11) < 1e-9 && std::abs(std::abs(cls.a21) - 1.0) < 1e-9, "SPICE CPW circuit element");
        const std::array<double, 1> wf{{10.0e9}};
        const auto wmode = rf::rectangular_waveguide_te10(22.86e-3, 10.16e-3, 1.0, 1.0, wf[0]);
        circuit::Circuit wgline;
        ta = wgline.node("ta");
        tb = wgline.node("tb");
        wgline.add_rectangular_waveguide_te10_line("WG", {ta, 0}, {tb, 0}, 22.86e-3, 10.16e-3, 1.0, 1.0, 0.03, wf,
                                                   wmode.wave_impedance);
        auto wgs = circuit::two_port_s_parameters(wgline, {"ta", "0"}, {"tb", "0"}, wf[0], wmode.wave_impedance);
        require(std::abs(wgs.a11) < 1e-9 && std::abs(std::abs(wgs.a21) - 1.0) < 1e-9,
                "SPICE rectangular-waveguide TE10 circuit element");
        circuit::Circuit sweep;
        in = sweep.node("in");
        out = sweep.node("out");
        sweep.add_voltage_source("V1", in, 0, 0.0, {1.0, 0.0});
        sweep.add_resistor("R1", in, out, 1000);
        sweep.add_resistor("R2", out, 0, 1000);
        auto dc = sweep.dc_sweep_voltage_source("V1", 0, 10, 3);
        require(dc.size() == 3 && near(sweep.voltage(dc.back().operating_point, "out"), 5.0, 1e-8),
                "SPICE DC source sweep");
        auto log = sweep.ac_log_sweep(10, 1e4, 4);
        require(log.size() == 4 && near(log.front().frequency_hz, 10) && near(log.back().frequency_hz, 1e4),
                "SPICE logarithmic AC sweep");
        circuit::Circuit noisy;
        out = noisy.node("out");
        noisy.add_resistor("R", out, 0, 1000);
        auto noise = noisy.output_noise(1e3, "out", 300);
        const double expected = std::sqrt(4.0 * 1.380649e-23 * 300.0 * 1000.0);
        require(std::abs(noise.output_noise_v_per_sqrt_hz / expected - 1.0) < 1e-6, "SPICE resistor thermal noise");
        const auto spectrum = circuit::output_noise_spectrum(noisy, "out", 10.0, 1.0e4, 9, 300.0);
        const double expected_integrated = expected * std::sqrt(1.0e4 - 10.0);
        require(spectrum.points.size() == 9 &&
                    std::abs(spectrum.integrated_rms_voltage / expected_integrated - 1.0) < 1e-6,
                "SPICE integrated noise sweep");
        circuit::Circuit psweep;
        auto psi = psweep.node("in"), pso = psweep.node("out");
        psweep.add_voltage_source("VP", psi, 0, 10.0);
        psweep.add_resistor("R1", psi, pso, 1000.0);
        psweep.add_resistor("R2", pso, 0, 1000.0);
        const std::array<double, 3> resistance_values{{500.0, 1000.0, 2000.0}};
        auto param = circuit::parameter_sweep_dc(
            psweep, resistance_values, [](circuit::Circuit& cc, double value) { cc.set_resistance("R2", value); });
        require(param.size() == 3 && near(psweep.voltage(param[0].operating_point, "out"), 10.0 / 3.0, 1e-8) &&
                    near(psweep.voltage(param[2].operating_point, "out"), 20.0 / 3.0, 1e-8),
                "SPICE general parameter sweep");
        circuit::Circuit thermal;
        auto tv = thermal.node("v"), td = thermal.node("d");
        thermal.add_voltage_source("VT", tv, 0, 1.0);
        thermal.add_resistor("RT", tv, td, 1000.0);
        thermal.add_diode_model("DT", {1e-14, 1.0, 300.0});
        thermal.add_diode("D", td, 0, "DT");
        auto tsweep = thermal.dc_temperature_sweep(280.0, 340.0, 2);
        require(tsweep.size() == 2 && tsweep[0].operating_point.converged && tsweep[1].operating_point.converged &&
                    std::abs(thermal.voltage(tsweep[0].operating_point, "d") -
                             thermal.voltage(tsweep[1].operating_point, "d")) > 1e-3,
                "SPICE device temperature sweep");
    }
    {
        const char* advanced = "VCTRL ctrl 0 1\nE1 out 0 ctrl 0 2\nR1 out 0 1k\n.model QN NPN (IS=1f BF=100 "
                               "BR=1)\n.model JN NJF (VTO=-2 BETA=1m)\n.model SM SW (RON=1 ROFF=1T VT=2.5 VH=.1)\n";
        auto parsed = circuit::Circuit::parse_spice(advanced);
        auto op = parsed.dc_operating_point();
        require(op.converged && near(parsed.voltage(op, "out"), 2.0, 1e-8), "SPICE controlled-source/model parser");
        const char* parametric = ".param RBASE=1k SCALE=2\nV1 in 0 10\nR1 in out {RBASE*SCALE}\nR2 out 0 RBASE\n";
        auto pc = circuit::Circuit::parse_spice(parametric);
        auto pop = pc.dc_operating_point();
        require(pop.converged && near(pc.voltage(pop, "out"), 10.0 / 3.0, 1e-8),
                "SPICE PARAM/expression netlist values");
        const std::unordered_map<std::string, double> ep{{"A", 2.0}, {"B", 3.0}};
        require(near(circuit::evaluate_spice_expression("{A^2+B*2}", ep), 10.0, 1e-12), "SPICE expression evaluator");
        const char* hierarchical = ".param RBASE=1k\n.subckt DIV IN OUT PARAMS: RUP=1k RDN=1k\nRUPPER IN MID "
                                   "{RUP}\nRLOWER MID OUT {RDN}\n.ends DIV\n.subckt WRAP A B PARAMS: SCALE=2\nXINNER A "
                                   "B DIV RUP=RBASE RDN={RBASE*SCALE}\n.ends WRAP\nVTOP in 0 9\nXW in 0 WRAP SCALE=3\n";
        auto hc = circuit::Circuit::parse_spice(hierarchical);
        auto hop = hc.dc_operating_point();
        require(hop.converged && near(hc.voltage(hop, "XW:XINNER:MID"), 6.75, 1e-8),
                "SPICE nested SUBCKT parameter passing");
    }
    {
        circuit::Circuit sense;
        auto in = sense.node("in"), out = sense.node("out");
        sense.add_voltage_source("V1", in, 0, 10.0);
        sense.add_resistor("R1", in, out, 1000);
        sense.add_resistor("R2", out, 0, 1000);
        auto sensitivity = circuit::voltage_source_dc_sensitivity(sense, "V1", "out", 1e-4);
        require(near(sensitivity.derivative, 0.5, 1e-8), "SPICE numerical DC sensitivity");
        const std::array<circuit::ResistorTolerance, 2> tol{{{"R1", 0.0}, {"R2", 0.0}}};
        auto mc = circuit::monte_carlo_dc_resistors(sense, "out", tol, 8, 7);
        require(near(mc.mean, 5.0, 1e-8) && mc.standard_deviation < 1e-12, "SPICE Monte Carlo tolerancing baseline");
        std::vector<circuit::TransientPoint> waveform;
        constexpr std::size_t samples = 201;
        for (std::size_t i = 0; i < samples; ++i) {
            const double t = static_cast<double>(i) / 200.0;
            circuit::TransientPoint p;
            p.time = t;
            p.node_voltage = {0.0, std::sin(2.0 * std::numbers::pi * 10.0 * t) +
                                       0.1 * std::sin(2.0 * std::numbers::pi * 20.0 * t)};
            waveform.push_back(std::move(p));
        }
        auto f = circuit::fourier_measurement(waveform, 1, 10.0);
        require(std::abs(f.magnitude - 1.0) < 1e-10, "SPICE Fourier measurement");
        require(std::abs(circuit::total_harmonic_distortion(waveform, 1, 10.0, 3) - 0.1) < 1e-10,
                "SPICE THD measurement");
        const double vt = 1.380649e-23 * 300 / 1.602176634e-19;
        const double limited = circuit::limit_pn_junction_voltage(2.0, 0.6, vt, 1e-14);
        require(limited < 2.0 && limited > 0.6, "SPICE PN junction voltage limiting");
        circuit::Circuit behavioral;
        auto x = behavioral.node("x");
        behavioral.add_current_source("I", 0, x, 1e-3, {1e-3, 0});
        behavioral.add_static_device("XNL", {x}, [](std::span<const double> v) {
            const double current = 1e-3 * (v[0] + v[0] * v[0] * v[0]);
            const double g = 1e-3 * (1.0 + 3.0 * v[0] * v[0]);
            return circuit::StaticDeviceEvaluation{{current}, {g}};
        });
        auto bop = behavioral.dc_operating_point_homotopy();
        require(bop.converged && std::abs(bop.node_voltage[x] +
                                          bop.node_voltage[x] * bop.node_voltage[x] * bop.node_voltage[x] - 1.0) < 1e-8,
                "SPICE generic compact-device DC residual/Jacobian");
        auto bac = behavioral.ac(1e3, &bop);
        const double expected_ac = 1.0 / (1.0 + 3.0 * bop.node_voltage[x] * bop.node_voltage[x]);
        require(std::abs(std::abs(bac.node_voltage[x]) - expected_ac) < 1e-8,
                "SPICE generic compact-device AC Jacobian");
        circuit::Circuit oneport;
        auto port_node = oneport.node("p");
        oneport.add_resistor("RLOAD", port_node, 0, 100.0);
        const std::array<circuit::CircuitPort, 1> ports{{{port_node, 0, 50.0}}};
        auto extracted = circuit::extract_s_parameters(oneport, ports, 1e9);
        require(std::abs(extracted(0, 0) - rf::Complex{1.0 / 3.0, 0.0}) < 1e-10,
                "SPICE circuit one-port S-parameter extraction");
    }
}

} // namespace
int main() {
    try {
        test_runtime();
        test_fvm_field();
        test_projection();
        test_io();
        test_optics();
        test_electrochem();
        test_extended_physics();
        test_workflow_and_diagnostics();
        test_hpc_extensions();
        test_fdtd_chemistry_extensions();
        test_analysis_chemistry_workflow();
        test_cfd_framework_extensions();
        test_fem_fdtd_material_extensions();
        test_coupled_cfd_extensions();
        test_lbm_advanced_extensions();
        test_cpu_solver_gap_closures();
        test_next_phase_portable();
        test_portable_physics_depth();
        test_rf_and_spice();
        std::cout << "completion tests passed\n";
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "completion test failed: " << e.what() << '\n';
        return 1;
    }
}
