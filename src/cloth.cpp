#include <iostream>
#include <math.h>
#include <random>
#include <vector>

#include "cloth.h"
#include "collision/plane.h"
#include "collision/sphere.h"

using namespace std;

Cloth::Cloth(double width, double height, int num_width_points,
             int num_height_points, float thickness) {
  this->width = width;
  this->height = height;
  this->num_width_points = num_width_points;
  this->num_height_points = num_height_points;
  this->thickness = thickness;

  buildGrid();
  buildClothMesh();
}

Cloth::~Cloth() {
  point_masses.clear();
  springs.clear();

  if (clothMesh) {
    delete clothMesh;
  }
}

void Cloth::buildGrid() {
  // TODO (Part 1): Build a grid of masses and springs.
    double w_increment = width / num_width_points;
    double h_increment = height / num_height_points;
    for (int j = 0; j < num_height_points; j++) {
        for (int i = 0; i < num_width_points; i++) {
            Vector3D pos;
            if (orientation == 0) {
                pos = Vector3D(w_increment * i, 1.0, h_increment * j);
            } else {
                double offset = (double)rand() / (double)RAND_MAX;
                double zVal = (-1.0/1000.0) + offset * (2.0/1000.0);
                cout << zVal << endl;
                pos = Vector3D(w_increment * i, h_increment * j, zVal);
            }
            std::vector<int> coord = {i,j};
            bool pinVal = std::count(pinned.begin(), pinned.end(), coord) > 0;
            PointMass pm = PointMass(pos, pinVal);
            point_masses.push_back(pm);
        }
    }
    for (int i = 0; i < num_height_points; i++) {
        for (int j = 0; j < num_width_points; j++) {
            int index = i * num_width_points + j;
            if (j > 0) {
                springs.push_back(Spring(&point_masses[index], &point_masses[index-1], STRUCTURAL));
            }
            if (i > 0) {
                springs.push_back(Spring(&point_masses[index], &point_masses[index-num_width_points], STRUCTURAL));
            }
            if (j > 1) {
                springs.push_back(Spring(&point_masses[index], &point_masses[index-2], BENDING));
            }
            if (i > 1) {
                springs.push_back(Spring(&point_masses[index], &point_masses[index-2*num_width_points], BENDING));
            }
            if (i > 0 && j > 0) {
                springs.push_back(Spring(&point_masses[index], &point_masses[index-num_width_points-1], SHEARING));
            }
            if (i > 0 && j < num_width_points - 1) {
                springs.push_back(Spring(&point_masses[index], &point_masses[index-num_width_points+1], SHEARING));
            }
        }
    }
    
}

void Cloth::simulate(double frames_per_sec, double simulation_steps, ClothParameters *cp,
                     vector<Vector3D> external_accelerations,
                     vector<CollisionObject *> *collision_objects) {
  double mass = width * height * cp->density / num_width_points / num_height_points;
  double delta_t = 1.0f / frames_per_sec / simulation_steps;
    
    

    
  // TODO (Part 2): Compute total force acting on each point mass.
    Vector3D netAccel = Vector3D();
    for (auto v : external_accelerations) {
        netAccel += v;
    }
    Vector3D netForce = mass * netAccel;
    
    for (PointMass &p : point_masses) {
        p.forces = netForce;
    }
    for (Spring &s : springs) {
        Vector3D curSpring = s.pm_b->position - s.pm_a->position;
        double magnitude = curSpring.norm();
        double springForce = cp->ks * (magnitude - s.rest_length);
        curSpring.normalize();
        s.pm_a->forces += springForce * curSpring;
        s.pm_b->forces -= springForce * curSpring;
    }

  // TODO (Part 2): Use Verlet integration to compute new point mass positions
    for (PointMass &p : point_masses) {
        if (!p.pinned) {
            Vector3D newPosition = p.position + (1.0 - cp->damping / 100.0) * (p.position - p.last_position) + p.forces / mass * pow(delta_t, 2);
            p.last_position = p.position;
            p.position = newPosition;
        }
    }
    
  // TODO (Part 4): Handle self-collisions.
    build_spatial_map();
    for (PointMass &p : point_masses) {
        self_collide(p, simulation_steps);
    }
    


  // TODO (Part 3): Handle collisions with other primitives.
    for (PointMass &p : point_masses) {
        for (CollisionObject *c : *collision_objects) {
            c->collide(p);
        }
    }


  // TODO (Part 2): Constrain the changes to be such that the spring does not change
  // in length more than 10% per timestep [Provot 1995].
    for (Spring &s : springs) {
        Vector3D curSpring = s.pm_b->position - s.pm_a->position;
        if (curSpring.norm() - s.rest_length * 1.1 > 0) {
            double offset = curSpring.norm() - s.rest_length * 1.1;
            Vector3D a_offset = Vector3D();
            Vector3D b_offset = Vector3D();
            if (!s.pm_a->pinned && !s.pm_b->pinned) {
                a_offset = curSpring / curSpring.norm() * offset * 0.5;
                b_offset = -curSpring / curSpring.norm() * offset * 0.5;
            } else if (!s.pm_b->pinned && s.pm_a->pinned) {
                b_offset = -curSpring / curSpring.norm() * offset;
            } else if (!s.pm_a->pinned && s.pm_b->pinned) {
                a_offset = curSpring / curSpring.norm() * offset;
            }
            s.pm_a->position += a_offset;
            s.pm_b->position += b_offset;
        }
    }
}

void Cloth::build_spatial_map() {
  for (const auto &entry : map) {
    delete(entry.second);
  }
  map.clear();

  // TODO (Part 4): Build a spatial map out of all of the point masses.
    for (PointMass &p : point_masses) {
        float hash = hash_position(p.position);
        if (map.count(hash) == 0) {
            map[hash] = new std::vector<PointMass *>();
        }
        map[hash]->push_back(&p);
    }
}

void Cloth::self_collide(PointMass &pm, double simulation_steps) {
  // TODO (Part 4): Handle self-collision for a given point mass.
    Vector3D total = Vector3D();
    int count = 0;
    float hash = hash_position(pm.position);
    if (map.count(hash)) {
        for (PointMass *p : *map[hash]) {
            if ((p->position - pm.position).norm() != 0) {
                double l = (p->position - pm.position).norm();
                if (l < 2.0 * thickness) {
                    count++;
                    double offset = 2.0 * thickness - l;
                    Vector3D direction = pm.position - p->position;
                    direction.normalize();
                    total += direction * offset;
                }
            }
        }
    }
    if (count) {
        total = total / (double)count / simulation_steps;
        pm.position += total;
    }
}

float Cloth::hash_position(Vector3D pos) {
  // TODO (Part 4): Hash a 3D position into a unique float identifier that represents membership in some 3D box volume.
    double w = 3.0 * width / num_width_points;
    double h = 3.0 * height / num_height_points;
    double t = max(w, h);
    double x_partition = pos.x - fmod(pos.x, w);
    double y_partition = pos.y - fmod(pos.y, h);
    double z_partition = pos.z - fmod(pos.z, t);
    return (float)(pow((x_partition), 1) + pow((y_partition), 2) + pow((z_partition), 3));

}

///////////////////////////////////////////////////////
/// YOU DO NOT NEED TO REFER TO ANY CODE BELOW THIS ///
///////////////////////////////////////////////////////

void Cloth::reset() {
  PointMass *pm = &point_masses[0];
  for (int i = 0; i < point_masses.size(); i++) {
    pm->position = pm->start_position;
    pm->last_position = pm->start_position;
    pm++;
  }
}

void Cloth::buildClothMesh() {
  if (point_masses.size() == 0) return;

  ClothMesh *clothMesh = new ClothMesh();
  vector<Triangle *> triangles;

  // Create vector of triangles
  for (int y = 0; y < num_height_points - 1; y++) {
    for (int x = 0; x < num_width_points - 1; x++) {
      PointMass *pm = &point_masses[y * num_width_points + x];
      // Get neighboring point masses:
      /*                      *
       * pm_A -------- pm_B   *
       *             /        *
       *  |         /   |     *
       *  |        /    |     *
       *  |       /     |     *
       *  |      /      |     *
       *  |     /       |     *
       *  |    /        |     *
       *      /               *
       * pm_C -------- pm_D   *
       *                      *
       */
      
      float u_min = x;
      u_min /= num_width_points - 1;
      float u_max = x + 1;
      u_max /= num_width_points - 1;
      float v_min = y;
      v_min /= num_height_points - 1;
      float v_max = y + 1;
      v_max /= num_height_points - 1;
      
      PointMass *pm_A = pm                       ;
      PointMass *pm_B = pm                    + 1;
      PointMass *pm_C = pm + num_width_points    ;
      PointMass *pm_D = pm + num_width_points + 1;
      
      Vector3D uv_A = Vector3D(u_min, v_min, 0);
      Vector3D uv_B = Vector3D(u_max, v_min, 0);
      Vector3D uv_C = Vector3D(u_min, v_max, 0);
      Vector3D uv_D = Vector3D(u_max, v_max, 0);
      
      
      // Both triangles defined by vertices in counter-clockwise orientation
      triangles.push_back(new Triangle(pm_A, pm_C, pm_B, 
                                       uv_A, uv_C, uv_B));
      triangles.push_back(new Triangle(pm_B, pm_C, pm_D, 
                                       uv_B, uv_C, uv_D));
    }
  }

  // For each triangle in row-order, create 3 edges and 3 internal halfedges
  for (int i = 0; i < triangles.size(); i++) {
    Triangle *t = triangles[i];

    // Allocate new halfedges on heap
    Halfedge *h1 = new Halfedge();
    Halfedge *h2 = new Halfedge();
    Halfedge *h3 = new Halfedge();

    // Allocate new edges on heap
    Edge *e1 = new Edge();
    Edge *e2 = new Edge();
    Edge *e3 = new Edge();

    // Assign a halfedge pointer to the triangle
    t->halfedge = h1;

    // Assign halfedge pointers to point masses
    t->pm1->halfedge = h1;
    t->pm2->halfedge = h2;
    t->pm3->halfedge = h3;

    // Update all halfedge pointers
    h1->edge = e1;
    h1->next = h2;
    h1->pm = t->pm1;
    h1->triangle = t;

    h2->edge = e2;
    h2->next = h3;
    h2->pm = t->pm2;
    h2->triangle = t;

    h3->edge = e3;
    h3->next = h1;
    h3->pm = t->pm3;
    h3->triangle = t;
  }

  // Go back through the cloth mesh and link triangles together using halfedge
  // twin pointers

  // Convenient variables for math
  int num_height_tris = (num_height_points - 1) * 2;
  int num_width_tris = (num_width_points - 1) * 2;

  bool topLeft = true;
  for (int i = 0; i < triangles.size(); i++) {
    Triangle *t = triangles[i];

    if (topLeft) {
      // Get left triangle, if it exists
      if (i % num_width_tris != 0) { // Not a left-most triangle
        Triangle *temp = triangles[i - 1];
        t->pm1->halfedge->twin = temp->pm3->halfedge;
      } else {
        t->pm1->halfedge->twin = nullptr;
      }

      // Get triangle above, if it exists
      if (i >= num_width_tris) { // Not a top-most triangle
        Triangle *temp = triangles[i - num_width_tris + 1];
        t->pm3->halfedge->twin = temp->pm2->halfedge;
      } else {
        t->pm3->halfedge->twin = nullptr;
      }

      // Get triangle to bottom right; guaranteed to exist
      Triangle *temp = triangles[i + 1];
      t->pm2->halfedge->twin = temp->pm1->halfedge;
    } else {
      // Get right triangle, if it exists
      if (i % num_width_tris != num_width_tris - 1) { // Not a right-most triangle
        Triangle *temp = triangles[i + 1];
        t->pm3->halfedge->twin = temp->pm1->halfedge;
      } else {
        t->pm3->halfedge->twin = nullptr;
      }

      // Get triangle below, if it exists
      if (i + num_width_tris - 1 < 1.0f * num_width_tris * num_height_tris / 2.0f) { // Not a bottom-most triangle
        Triangle *temp = triangles[i + num_width_tris - 1];
        t->pm2->halfedge->twin = temp->pm3->halfedge;
      } else {
        t->pm2->halfedge->twin = nullptr;
      }

      // Get triangle to top left; guaranteed to exist
      Triangle *temp = triangles[i - 1];
      t->pm1->halfedge->twin = temp->pm2->halfedge;
    }

    topLeft = !topLeft;
  }

  clothMesh->triangles = triangles;
  this->clothMesh = clothMesh;
}
