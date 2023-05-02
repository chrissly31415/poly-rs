// visualize.rs
use crate::{Array3, Molecule, MAX_TYPE};

use kiss3d::light::Light;
use kiss3d::nalgebra::{Point2, Point3, Point4, Translation3, UnitQuaternion, Vector3};
use kiss3d::scene::SceneNode;
use kiss3d::text::Font;
use kiss3d::window::Window;
use std::cell::RefCell;
use std::collections::{HashMap, HashSet};
use std::path::Path;
use std::rc::Rc;

use petgraph::graph::Graph;
use petgraph::Undirected;

fn draw_axes_and_labels(window: &mut Window) {
    //let origin = Point3::<f32>::origin();
    let axis_length = 1.0;
    let axes_rot = [Vector3::z_axis(), Vector3::y_axis(), Vector3::x_axis()];
    let axes_shift = [Vector3::x_axis(), Vector3::y_axis(), Vector3::z_axis()];

    let axes_colors = [
        Point3::new(1.0, 0.0, 0.0), // X Red
        Point3::new(0.0, 1.0, 0.0), // Y Green is default, no rotation
        Point3::new(0.0, 0.0, 1.0), // Z Blue
    ];

    for i in 0..3 {
        let mut axis = window.add_cylinder(0.05, axis_length);
        axis.set_color(axes_colors[i].x, axes_colors[i].y, axes_colors[i].z);
        if i != 1 {
            let rotation = UnitQuaternion::from_axis_angle(&axes_rot[i], 1.5707964);
            axis.set_local_rotation(rotation)
        };
        axis.set_local_translation(Translation3::from(axes_shift[i].into_inner() * 0.5));
    }

    //label
    /*     let label_offset = 0.2;
    let x_label_position = Point2::new(axes_rot[0].x + label_offset, axes_rot[0].y);
    let y_label_position = Point2::new(axes_rot[1].x, axes_rot[1].y + label_offset);
    let z_label_position = Point2::new(axes_rot[2].x, axes_rot[2].z + label_offset);

    let label_color = Point3::new(1.0, 1.0, 1.0); // White
    let label_scale = 0.05;

    let font_path = Path::new("/home/loschen/calc/poly-rs/src/IBM_Plex_Sans/IBMPlexSans-Light.ttf"); // Update the path to your font file
    let font = Font::new(font_path).unwrap_or_else(|| {
        panic!("Failed to load font from file");
    });

    window.draw_text("X", &x_label_position, label_scale, &font, &label_color);
    window.draw_text("Y", &y_label_position, label_scale, &font, &label_color);
    window.draw_text("Z", &z_label_position, label_scale, &font, &label_color); */
}

use petgraph::algo::kosaraju_scc;
pub fn visualize_chains(graph: &Graph<Rc<RefCell<Molecule>>, (), Undirected>) {
    let mut window: Window = Window::new("Chain Viewer");
    window.set_light(Light::StickToCamera);

    draw_axes_and_labels(&mut window);

    let scc = kosaraju_scc(graph);

    let mut color_map = HashMap::new();
    let mut unique_lengths: HashSet<usize> = HashSet::new();

    for chain in &scc {
        unique_lengths.insert(chain.len());
    }

    // Sort unique_lengths in descending order
    let mut unique_lengths_sorted: Vec<usize> = unique_lengths.into_iter().collect();
    unique_lengths_sorted.sort_unstable_by(|a, b| b.cmp(a));

    let num_colors = unique_lengths_sorted.len();
    let mut idx = 0;

    for len in unique_lengths_sorted.iter() {
        let r = (idx * 255 / num_colors) as f32 / 255.0;
        let g = ((num_colors - idx) * 255 / num_colors) as f32 / 255.0;
        let b = 0.0;
        let alpha = 0.5;
        color_map.insert(*len, Point4::new(r, g, b, alpha));
        idx += 1;
    }

    for chain in scc.iter() {
        let chain_color = color_map.get(&chain.len()).unwrap();
        for node_idx in chain {
            let molecule = graph.node_weight(*node_idx).unwrap();
            let molecule_borrowed = molecule.borrow();
            let (i, j, k) = molecule_borrowed.position;

            let mut node = window.add_sphere(1.0);
            node.set_local_translation(Translation3::from(Vector3::new(
                i as f32, j as f32, k as f32,
            )));
            node.set_color(chain_color.x, chain_color.y, chain_color.z);
        }
    }

    while window.render() {}
}

pub fn visualize_layers(
    grid: &Array3<Option<Rc<RefCell<Molecule>>>>,
    molecules: &Vec<Rc<RefCell<Molecule>>>,
    rad: f32,
    show_bonds: bool,
    use_bfm: bool,
    show_grid: bool,
) {
    let mut window: Window = Window::new("Molecule Viewer");
    window.set_light(Light::StickToCamera);
    let alpha = 0.5;

    draw_axes_and_labels(&mut window);

    let color_map: Vec<Point4<f32>> = (0..MAX_TYPE)
        .map(|i| {
            let r = (i * 255 / MAX_TYPE) as f32 / 255.0;
            let g = ((MAX_TYPE - i) * 255 / MAX_TYPE) as f32 / 255.0;
            let b = 0.0;
            Point4::new(r, g, b, alpha) // Add the alpha component to the color
        })
        .collect();

    if show_grid {
        for ((i, j, k), _) in grid.indexed_iter() {
            let mut grid_node = window.add_sphere(0.1); // Small white sphere for grid points
            grid_node.set_color(1.0, 1.0, 1.0);
            grid_node.set_local_translation(Translation3::from(Vector3::new(
                i as f32, j as f32, k as f32,
            )));
        }
    }

    let cube_size = 1.0;
    let nodes: Vec<SceneNode> = molecules
        .iter()
        .map(|molecule| {
            let molecule_borrowed = molecule.borrow();
            let (i, j, k) = molecule_borrowed.position;
            let color = color_map[molecule_borrowed.functional_group.id as usize];

            let mut node;
            if use_bfm {
                node = window.add_cube(cube_size, cube_size, cube_size);
                node.set_local_translation(Translation3::from(Vector3::new(
                    (i as f32) + cube_size / 2.0,
                    (j as f32) + cube_size / 2.0,
                    (k as f32) + cube_size / 2.0,
                )));
            } else {
                node = window.add_sphere(rad);
                node.set_local_translation(Translation3::from(Vector3::new(
                    i as f32, j as f32, k as f32,
                )));
            }
            node.set_color(color.x, color.y, color.z);

            //node.set_alpha(color.w); // Set the alpha component for transparency

            node
        })
        .collect();

    if show_bonds {
        // Create bond cylinders
        let mut bond_cylinders: Vec<SceneNode> = vec![];
        let mut drawn_bonds = std::collections::HashSet::new();

        for molecule in molecules {
            let molecule_borrowed = molecule.borrow();
            for neighbor in &molecule_borrowed.neighbors {
                let neighbor_borrowed = neighbor.borrow();
                let bond_key = if Rc::ptr_eq(molecule, neighbor) {
                    (neighbor.as_ptr(), molecule.as_ptr())
                } else {
                    (molecule.as_ptr(), neighbor.as_ptr())
                };

                if !drawn_bonds.contains(&bond_key) {
                    drawn_bonds.insert(bond_key);

                    let (i, j, k) = molecule_borrowed.position;
                    let (ni, nj, nk) = neighbor_borrowed.position;

                    let p1;
                    let p2;
                    if use_bfm {
                        p1 = Vector3::new(
                            (i as f32) + cube_size / 2.0,
                            (j as f32) + cube_size / 2.0,
                            (k as f32) + cube_size / 2.0,
                        );
                        p2 = Vector3::new(
                            (ni as f32) + cube_size / 2.0,
                            (nj as f32) + cube_size / 2.0,
                            (nk as f32) + cube_size / 2.0,
                        );
                    } else {
                        p1 = Vector3::new(i as f32, j as f32, k as f32);
                        p2 = Vector3::new(ni as f32, nj as f32, nk as f32);
                    }

                    let bond_length = (p2 - p1).norm();

                    let bond_direction = (p2 - p1).normalize();
                    let rotation =
                        UnitQuaternion::rotation_between(&Vector3::z_axis(), &bond_direction)
                            .unwrap_or(UnitQuaternion::identity());
                    let midpoint = (p1 + p2) / 2.0;
                    let mut cylinder = window.add_cylinder(1.0, 1.0);
                    cylinder.set_local_translation(Translation3::from(midpoint));
                    cylinder.set_color(1.0, 0.0, 0.0);
                    cylinder.set_local_scale(0.1, 0.1, bond_length);

                    cylinder.set_local_rotation(rotation);
                    bond_cylinders.push(cylinder);
                }
            }
        }
    }

    // Render loop
    while window.render() {}
    window.close();
}
