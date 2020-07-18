#include "driver_state.h"
#include <cstring>
#include <algorithm>
#include <climits>
#include <cfloat>
#include <vector>


#define C_MAX 255
#define OUT 1
#define IN 2

driver_state::driver_state()
{
}

driver_state::~driver_state()
{
    delete [] image_color;
    delete [] image_depth;
}

// This function should allocate and initialize the arrays that store color and
// depth.  This is not done during the constructor since the width and height
// are not known when this class is constructed.
void initialize_render(driver_state& state, int width, int height)
{
    state.image_width=width;
    state.image_height=height;
    state.image_color=0;
    state.image_depth=0;
    std::cout<<"TODO: allocate and initialize state.image_color and state.image_depth."<<std::endl;

    int len = width * height;

    state.image_color = new pixel[len];
    state.image_depth = new float[len];

    for (int i = 0; i < len; i++) {
        state.image_color[i] = make_pixel(0, 0, 0);
        state.image_depth[i] = FLT_MAX;
    }
}

// This function will be called to render the data that has been stored in this class.
// Valid values of type are:
//   render_type::triangle - Each group of three vertices corresponds to a triangle.
//   render_type::indexed -  Each group of three indices in index_data corresponds
//                           to a triangle.  These numbers are indices into vertex_data.
//   render_type::fan -      The vertices are to be interpreted as a triangle fan.
//   render_type::strip -    The vertices are to be interpreted as a triangle strip.
void render(driver_state& state, render_type type)
{
    // std::cout<<"TODO: implement rendering."<<std::endl;
    int triangles;
    int vert_index = 0;
    
    data_geometry * data_geos = new data_geometry[3];

    switch (type) {
    case render_type::triangle:
        triangles = state.num_vertices / 3;
        
        for (int i = 0; i < triangles; i++) {
            for (int i = 0; i < 3; i++) {
                data_geos[i].data = state.vertex_data + vert_index;
                vert_index += state.floats_per_vertex;
            }
            data_vertex data_vert;
            for (int i = 0; i < 3; i++) {
                data_vert.data = data_geos[i].data;
                state.vertex_shader(data_vert, data_geos[i], state.uniform_data);
            }
            clip_triangle(state, (const data_geometry **)(&data_geos), 0);
        }
        break;

    case render_type::indexed:
        for (int i = 0; i < state.num_triangles; i++) {
            for (int i = 0; i < 3; i++) {
                data_geos[i].data = state.vertex_data + state.index_data[vert_index] * state.floats_per_vertex;
                vert_index++;
            }
            data_vertex data_vert;
            for (int i = 0; i < 3; i++) {
                data_vert.data = data_geos[i].data;
                state.vertex_shader(data_vert, data_geos[i], state.uniform_data);
            }
            clip_triangle(state, (const data_geometry **)(&data_geos), 0);
        }
        break;

    case render_type::fan:
        triangles = state.num_vertices - 2;
        vert_index = 1;
        data_geos[0].data = state.vertex_data;
        for (int i = 0; i < triangles; i++) {
            for (int i = 1; i < 3; i++) {
                data_geos[i].data = state.vertex_data + (vert_index * state.floats_per_vertex);
                vert_index++;
            }
            vert_index--;
            data_vertex data_vert;
            for (int i = 0; i < 3; i++) {
                data_vert.data = data_geos[i].data;
                state.vertex_shader(data_vert, data_geos[i], state.uniform_data);
            }
            clip_triangle(state, (const data_geometry **)(&data_geos), 0);
        }       
        break;

    case render_type::strip:
        triangles = state.num_vertices - 2;
        for (int i = 0; i < 3; i++) {
            data_geos[i].data = state.vertex_data + vert_index;
            vert_index += state.floats_per_vertex;
        }
        
        data_vertex data_vert;
        for (int i = 0; i < 3; i++) {
            data_vert.data = data_geos[i].data;
            state.vertex_shader(data_vert, data_geos[i], state.uniform_data);
        }
        clip_triangle(state, (const data_geometry **)(&data_geos), 0);
        for (int i = 1; i < triangles; i++) {
            if (i % 2) {
                data_geos[0].data = data_geos[1].data;
                data_geos[1].data = state.vertex_data + vert_index;
            } else {
                data_geos[0].data = data_geos[2].data;
                data_geos[2].data = state.vertex_data + vert_index;
            }

            vert_index += state.floats_per_vertex;
            data_vertex data_vert;
            for (int i = 0; i < 3; i++) {
                data_vert.data = data_geos[i].data;
                state.vertex_shader(data_vert, data_geos[i], state.uniform_data);
            }
            clip_triangle(state, (const data_geometry **)(&data_geos), 0);
        }
        break;

    default:
        break;
    }

    delete[] data_geos;

}

// This function clips a triangle (defined by the three vertices in the "in" array).
// It will be called recursively, once for each clipping face (face=0, 1, ..., 5) to
// clip against each of the clipping faces in turn.  When face=6, clip_triangle should
// simply pass the call on to rasterize_triangle.
void clip_triangle(driver_state& state, const data_geometry* in[3],int face)
{
    if(face==6)
    {
        rasterize_triangle(state, in);
        return;
    }
    // std::cout<<"TODO: implement clipping. (The current code passes the triangle through without clipping them.)"<<std::endl;
    std::vector<data_geometry *> tris;
    int sign = 2 * (face % 2) - 1;
    unsigned axis = face % 3;
    bool inside[3] = {0};
    
    tris.push_back(new data_geometry[3]);
    for (unsigned i = 0; i < 3; i++) {
        tris[tris.size() - 1][i].data = new float[MAX_FLOATS_PER_VERTEX];
        for (unsigned j = 0; j < 4; j++) {
            tris[tris.size() - 1][i].gl_Position[j] = (*in)[i].gl_Position[j];
        }

        for (int j = 0; j < state.floats_per_vertex; j++) {
            tris[tris.size() - 1][i].data[j] = (*in)[i].data[j];
        }
    }

    for (unsigned i = 0; i < 3; i++) {
        if (sign > 0) {
            inside[i] = (*in)[i].gl_Position[axis] <= (*in)[i].gl_Position[3];
        }
        if (sign < 0) {
            inside[i] = (*in)[i].gl_Position[axis] >= -1 * (*in)[i].gl_Position[3];
        }
    }

    bool track = true;
    for (unsigned i = 0; i < 3; i++) {
        if (!inside[i]) {
            track = false;
        }
    }
    if (!track) {
        if (inside[0] && !inside[1] && !inside[2]) { 
            createT(tris, axis, sign, 0, 1, 2, state, OUT);
        } else if (!inside[0] && inside[1] && !inside[2]) { 
            createT(tris, axis, sign, 1, 2, 0, state, OUT);
        } else if (!inside[0] && !inside[1] && inside[2]) { 
            createT(tris, axis, sign, 2, 0, 1, state, OUT);
        } else if (!inside[0] && inside[1] && inside[2]) {
            createT(tris, axis, sign, 0, 1, 2, state, IN);
        } else if (inside[0] && !inside[1] && inside[2]) {
            createT(tris, axis, sign, 1, 2, 0, state, IN);
        } else if (inside[0] && inside[1] && !inside[2]) {
            createT(tris, axis, sign, 2, 0, 1, state, IN);
        }
    }

    for (unsigned i = 0; i < tris.size(); i++) {
        clip_triangle(state,(const data_geometry **)(&(tris[i])),face+1);
    }

}

void createT(std::vector<data_geometry *>& tris, unsigned axis, int sign, unsigned ind0, unsigned ind1, unsigned ind2, const driver_state& state, const int at) {
    data_geometry a = tris[0][ind0];
    data_geometry b = tris[0][ind1];
    data_geometry c = tris[0][ind2];

    float bNum = sign * b.gl_Position[3] - b.gl_Position[axis];
    float abWeight = bNum / ((a.gl_Position[axis] - sign * a.gl_Position[3]) + bNum);
    float cNum = sign * c.gl_Position[3] - c.gl_Position[axis];
    float acWeight = cNum / ((a.gl_Position[axis] - sign * a.gl_Position[3]) + cNum);

    float abN = abWeight * tris[0][ind0].gl_Position[3] * 1.0f / (abWeight * tris[0][ind0].gl_Position[3] + (1 - abWeight) * tris[0][ind1].gl_Position[3]);
    float acN = acWeight * tris[0][ind0].gl_Position[3] * 1.0f / (acWeight * tris[0][ind0].gl_Position[3] + (1 - acWeight) * tris[0][ind2].gl_Position[3]);

    vec4 ab = abWeight * a.gl_Position + (1 - abWeight) * b.gl_Position;
    vec4 ac = acWeight * a.gl_Position + (1 - acWeight) * c.gl_Position;

    data_geometry * geos = new data_geometry[3];
    data_geometry * geos2 = new data_geometry[3];

    for (unsigned j = 0; j < 3; j++) {
        geos[j].data = new float[MAX_FLOATS_PER_VERTEX];
        if(at==IN){
            geos2[j].data = new float[MAX_FLOATS_PER_VERTEX];
        }
    }
    if(at==OUT){
        geos[0].gl_Position = a.gl_Position;
        geos[1].gl_Position = ab;
        geos[2].gl_Position = ac;
        tris.push_back(geos);
    }else{
        geos[0].gl_Position = ac;
        geos[1].gl_Position = b.gl_Position;
        geos[2].gl_Position = c.gl_Position;
        geos2[0].gl_Position = ab;
        geos2[1].gl_Position = b.gl_Position;
        geos2[2].gl_Position = ac;
        tris.push_back(geos);
        tris.push_back(geos2);
    }
    
    for (int i = 0; i < state.floats_per_vertex; i++) {
        tris[1][0].data[i] = tris[0][ind0].data[i];
        tris[1][1].data[i] = tris[0][ind1].data[i];
        tris[1][2].data[i] = tris[0][ind2].data[i];
        if(at==IN){
            tris[2][0].data[i] = tris[0][ind0].data[i];
            tris[2][1].data[i] = tris[0][ind1].data[i];
            tris[2][2].data[i] = tris[0][ind2].data[i];
        }
    }

    if(at==OUT){
        for (int j = 0; j < state.floats_per_vertex; j++) {
            tris[1][0].data[j] = tris[0][ind0].data[j];
        }
    }
    
    for (int i = 0; i < state.floats_per_vertex; i++) {
        switch (state.interp_rules[i]) {
        case interp_type::flat:
            for (unsigned j = 1; j < 3; j++) {
                tris[1][j].data[i] = tris[0][ind0].data[i];
                if(at==IN){
                    tris[2][j].data[i] = tris[0][ind0].data[i];
                }
            }
            break;
    
        case interp_type::smooth:
            if(at==OUT){
                tris[1][1].data[i] = abWeight * tris[0][ind0].data[i] + (1 - abWeight) * tris[0][ind1].data[i];
                tris[1][2].data[i] = acWeight * tris[0][ind0].data[i] + (1 - acWeight) * tris[0][ind2].data[i];
            }else{
                tris[1][0].data[i] = acWeight * tris[0][ind0].data[i] + (1 - acWeight) * tris[0][ind2].data[i];
                tris[2][0].data[i] = abWeight * tris[0][ind0].data[i] + (1 - abWeight) * tris[0][ind1].data[i];
                tris[2][2].data[i] = tris[1][0].data[i];
            }
            break;

        case interp_type::noperspective:
            if(at==OUT){
                tris[1][1].data[i] = abN * tris[0][ind0].data[i] + (1 - abN) * tris[0][ind1].data[i];
                tris[1][2].data[i] = acN * tris[0][ind0].data[i] + (1 - acN) * tris[0][ind2].data[i];
            }else{
                tris[1][0].data[i] = acN * tris[0][ind0].data[i] + (1 - acN) * tris[0][ind2].data[i];
                tris[2][0].data[i] = abN * tris[0][ind0].data[i] + (1 - abN) * tris[0][ind1].data[i];
                tris[2][2].data[i] = tris[1][0].data[i];
            }
            break;

        default:
            break;
        }
    }

    for (unsigned i = 0; i < 3; i++) {
        delete[] tris[0][i].data;
    }
    delete[] tris[0];
    tris.erase(tris.begin() + 0);
}

// Rasterize the triangle defined by the three vertices in the "in" array.  This
// function is responsible for rasterization, interpolation of data to
// fragments, calling the fragment shader, and z-buffering.
void rasterize_triangle(driver_state& state, const data_geometry* in[3])
{
    unsigned pixel;

    float x[3];
    float y[3];

    float z[3];
    float depth;

    float min_x, min_y;
    float max_x, max_y;

    float pt1[3];
    float pt2[3];
    float pt3[3]; 
    float area;
    float bary[3];

    data_fragment frag;
    frag.data = new float[MAX_FLOATS_PER_VERTEX];
    
    for (int iter = 0; iter < 3; iter++) {
        float w2 = state.image_width / 2.0f;
        float h2 = state.image_height / 2.0f;
        
        x[iter] = (w2 * (*in)[iter].gl_Position[0] / (*in)[iter].gl_Position[3] + (w2 - .5f));
        y[iter] = (h2 * (*in)[iter].gl_Position[1] / (*in)[iter].gl_Position[3] + (h2 - .5f));
    }

    area = .5f * ((x[1] * y[2] - x[2] * y[1]) - (x[0] * y[2] - x[2] * y[0]) + (x[0] * y[1] - x[1] * y[0]));

    pt1[0] = x[1] * y[2] - x[2] * y[1];
    pt1[1] = x[2] * y[0] - x[0] * y[2];
    pt1[2] = x[0] * y[1] - x[1] * y[0];
    pt2[0] = y[1] - y[2];
    pt2[1] = y[2] - y[0];
    pt2[2] = y[0] - y[1];
    pt3[0] = x[2] - x[1];
    pt3[1] = x[0] - x[2];
    pt3[2] = x[1] - x[0];

    min_x = state.image_width - 1;
    min_y = state.image_height - 1;
    max_x = 0;
    max_y = 0;

    for (int i = 0; i < 3; i++) {
        min_x = std::min(min_x, x[i]);
        min_y = std::min(min_y, y[i]);
        max_x = std::max(max_x, x[i]);
        max_y = std::max(max_y, y[i]);
    }

    min_x = std::max(min_x, 0.0f);
    min_y = std::max(min_y, 0.0f);
    max_x = std::min(max_x, state.image_width - 1.0f);
    max_y = std::min(max_y, state.image_height - 1.0f);


    for (int y = min_y + 1; y < max_y + 1; y++) {
        for (int x = min_x + 1; x < max_x + 1; x++) {
            for (int vert = 0; vert < 3; vert++) {
                bary[vert] = .5f * (pt1[vert] + (pt2[vert] * x) + (pt3[vert] * y)) / area;
            }

            float count = 0;
            for (unsigned i = 0; i < 3; i++) {
                z[i] = (*in)[i].gl_Position[2] / (*in)[i].gl_Position[3];
                count += z[i] * bary[i];
            }
            depth = count;

            pixel = x + y * state.image_width;

            bool bWeight = true;
            for (int i = 0; i < 3; i++) {
                if (bary[i] < 0) {
                    bWeight = false;
                    break;
                }
            }
            if (bWeight && depth < state.image_depth[pixel]) {
                float wBary[3];
                data_output out;
                float k = 0;

                for (int i = 0; i < state.floats_per_vertex; i++) {
                    switch (state.interp_rules[i]) {
                    
                    case interp_type::flat:
                        frag.data[i] = (*in)[0].data[i];
                        break;

                    case interp_type::smooth:
                        k = 0;
                        for (unsigned j = 0; j < 3; j++) {
                            k += bary[j] / (*in)[j].gl_Position[3];
                        }
                        count =0;
                        for (unsigned j = 0; j < 3; j++) {
                            wBary[j] = bary[j] / ((*in)[j].gl_Position[3]* k);
                            count += wBary[j] * (*in)[j].data[i];
                        }
                        frag.data[i] = count;
                        break;

                    case interp_type::noperspective:
                        count =0;
                        for (unsigned j = 0; j < 3; j++) {
                            count += bary[j] * (*in)[j].data[i];
                        }
                        frag.data[i] = count;
                        break;

                    default:
                        break;
                    }
                }

                state.fragment_shader(frag, out, state.uniform_data);
                state.image_color[pixel] = make_pixel(out.output_color[0] * C_MAX, out.output_color[1] * C_MAX, out.output_color[2] * C_MAX);
                state.image_depth[pixel] = depth;
            }
        }
    }
}
