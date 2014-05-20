
#include <opencv2/opencv.hpp>
#include <queue>
#include "map_learning.h"
#include "utils/graphs/graph_utils.h"


// map_learning
///////////////////////////////////////////////////////////////////////////////

map_learning::map_learning() :
    stat(),
    regions(),
    rcenters()
{
}

map_learning::map_learning(const ConfigDictionary& cfg) :
    stat(),
    regions(),
    rcenters()
{
    cfg_init(cfg);
}

map_learning::~map_learning()
{
}

void map_learning::cfg_init(const ConfigDictionary& cfg)
{
    if (cfg.isDefined("source_layer_index")) {
        source_layer = cfg.getValueInt("source_layer_index", 0) - 1;  // backward compatibility!
    } else {
        cfg.getValue(source_layer, "source_layer", true);
    }
    cfg.getValue(stat_dim, "nb_size", true);

    // map update
    center_val_threshold = cfg.getValueDouble("center_val_threshold", 0.0);
    center_val_threshold_rel = cfg.getValueDouble("center_val_threshold_rel", 0.5);
    nb_val_threshold_rel = cfg.getValueDouble("nb_val_threshold_rel", 0.6);
    nbthresh_min = cfg.getValueDouble("nbthresh_min", 0.0);
    nbthresh_max = cfg.getValueDouble("nbthresh_max", 0.0);
    seq_min_intersection_percent = cfg.getValueDouble("seq_min_intersection_percent", 0.0);
    seq_max_intersection_percent = cfg.getValueDouble("seq_max_intersection_percent", 0.5);

    // finding maxima
    max_max = cfg.getValueInt("max_max", 4);
    max_val_threshold = cfg.getValueDouble("max_val_threshold", 0.01);
    if (!cfg.isDefined("max_val_threshold"))
        max_val_threshold = cfg.getValueDouble("min_update_count_percent", 0.01);
    individual_max = cfg.getValueBool("individual_max", false);
    max_sigma = cfg.getValueDouble("max_sigma", 0.0);
    max_nbhood_mask = cfg.getValueInt("max_nbhood_mask", 5);
    max_radius = cfg.getValueInt("max_radius", 2); 
}

void map_learning::reset(part_lib* library)
{
    library->get_regions(1, regions, rcenters);
    stat.clear();
}

void map_learning::dispose()
{
    stat.clear();
}

// Prepares res for update. Sets ALLOW_UPDATE_ATTR to all nodes.
void map_learning::prepare_for_update(layer1_result* res)
{
    if (res == nullptr) 
        return;

    res->set_attr(res->shape_nodes[source_layer].begin(),
        res->shape_nodes[source_layer].end(), ALLOW_UPDATE_ATTR);
}

void map_learning::update(layer1_result* res)
{
    typedef map<node*, set<int> > rmap_t;

    if (res == nullptr)
        return;

    if (regions.empty()) 
        throw new_libhop_exception("map_learning::update: Map learner not initialized with library");

    vector<node*>& s_nodes = res->shape_nodes[source_layer];
    rmap_t rmap;

    get_region_map(rmap, res, source_layer, regions, rcenters);
    for (vector<node*>::iterator iter = s_nodes.begin(); iter != s_nodes.end(); ++iter) {
        node* cn = *iter;

        do {
            layer1_data* cd = (layer1_data*)cn->data;
            ipoint2 cp(cd->x, cd->y);
        
            if (cn->is_attr_set(ALLOW_UPDATE_ATTR) && cd->vval() >= center_val_threshold) {
                double nbt_min, nbt_max;

                if (nbthresh_min >= nbthresh_max) {
                    nbt_min = cd->vval() * nb_val_threshold_rel;
                    nbt_max = cd->vval() * (2 - nb_val_threshold_rel);
                } else {
                    nbt_min = nbthresh_min;
                    nbt_max = nbthresh_max;
                }
                set<int>& cset = rmap[cn];
                int xoff = stat_dim/2 - cd->x;
                int yoff = stat_dim/2 - cd->y;
                vector<node*> nbhood;

                res->get_nodes_circular(nbhood, cn, stat_dim/2, stat_dim/2, ALLOW_UPDATE_ATTR);

                // Update map with neighbors
                for (vector<node*>::iterator niter = nbhood.begin(); niter != nbhood.end(); ++niter) {
                    node* nb = *niter;
                    layer1_data* nbd = (layer1_data*)nb->data;
                    double nbval = nbd->vval();
                    set<int>& nbset = rmap[nb];
                    int is = intersection_size(cset, nbset);
                    int mins = min<int>((int)nbset.size(), (int)cset.size());

                    if ((is >= seq_min_intersection_percent * mins && is <= seq_max_intersection_percent * mins) &&
                            (nbt_min <= nbval && nbval <= nbt_max) &&
                            cp.distance2(nbd->x, nbd->y) > sqr(max_radius + 0.5)) {
                        matrix<double>& im = stat[ipoint2(cd->m, nbd->m)];

                        if (im.empty()) im.resize(stat_dim, stat_dim, 0.0);
                        im(xoff + nbd->x, yoff + nbd->y) += nbval; // * ::pow(dist2, -0.125)
                    }
                }

            } // end of "if (center->is_attr_set(ALLOW_UPDATE_ATTR) && cdata->val >= center_val_threshold)"

            cn = cd->next;
            break; // Only the best!
        } while (cn != nullptr);
    }
        
}

// Save statistics to files
// file_template should contain two %d format characters
void map_learning::display_statistics(const string& file_template, double blur_sigma /* = 0.75 */) const
{
    char fname[1000];
    int blur_dim = 4*blur_sigma + 1;
    img* blur_mask = img::gaussian_mask(blur_dim, blur_dim, blur_sigma);

    for (stat_map_t::const_iterator iter = stat.begin(); iter != stat.end(); ++iter) {
        matrix<double> tmp(iter->second);
        img im;

        tmp.add_border(blur_dim/2, blur_dim/2, 0.0);
        convolve_matrix(im, tmp, *blur_mask);

        sprintf(fname, file_template.c_str(), iter->first.x, iter->first.y);
        im.save_jet_colormap(fname, -1000);
    }
    delete blur_mask;
}

// Get vector of maxima
// stat: matrix of statistics
// max_maxima: max number of maxima
// max_thresh: each max must be >= max_thresh; 
void map_learning::get_maxima(max_vector_t& maxima, const matrix<double>& stat, double max_thresh) const
{
    typedef pair<double, ipoint2> queue_item_t;
    typedef priority_queue<queue_item_t> queue_t;

    const int border = 5;

    static img blurmask = img::gaussian_mask2(border, border, 0.75);

    maxima.clear();
    matrix<double> m(stat);
    matrix<double> dist;
    queue_t queue;
    
    m.add_border(border, border, 0.0);
    m = convolve_matrix(m, blurmask);
    m.change_values_leq(max_thresh, 0.0);

    if (max_sigma > 0.0)
        dist = img::gaussian_mask2(max_nbhood_mask, max_nbhood_mask, max_sigma);

    for_each_xy_int (m, i, j) {
        double lmax = m.maximum_c((int)i, (int)j, max_radius, max_radius);

        if (lmax > 1E-6 && m(i, j) == lmax) 
            queue.push(queue_item_t(lmax, ipoint2(i, j)));
    }

    matrix<int> bm(m.width, m.height, 0);

    bm.set_region_c(m.width/2, m.height/2, max_radius, max_radius, 1);
    while (!queue.empty() && (int)maxima.size() < max_max) {
        ipoint2 p = queue.top().second;

        if (bm(p.x, p.y) == 0) {
            bm.set_region_c(p.x, p.y, max_radius, max_radius, 1);
            maxima.push_back(map_learning::max_item_t());
            maxima.back().pos = ipoint2(p.x - m.width/2, p.y - m.height/2);
            if (!dist.empty())
                maxima.back().dist = dist;
            else {
                int radius2 = max_nbhood_mask/2;

                get_region_sym(maxima.back().dist, m, p.x - radius2, p.x + radius2 + 1,
                    p.y - radius2, p.y + radius2 + 1);
                max_normalize(maxima.back().dist);
            }
        }
        queue.pop();
    }
}

void map_learning::get_maxima(max_map_t& result) const
{
    double max = 0.0;

    result.clear();
    if (!individual_max) {
        for (stat_map_t::const_iterator iter = stat.begin(); iter != stat.end(); ++iter) {
            const stat_t& m = iter->second;
            int r = m.maximum();

            if (r > max) max = r;
        }
        max *= max_val_threshold;
    }
    for (stat_map_t::const_iterator iter = stat.begin(); iter != stat.end(); ++iter) {
        const stat_t& m = iter->second;
        
        if (individual_max)
            max = max_val_threshold * m.maximum();
        get_maxima(result[iter->first], m, max);
    }
}

void map_learning::read_maxima_from_stream(max_map_t& maxima, vector<matrix<bool> >& regions, vector<ipoint2>& rcenters,
        const string& fname)
{
    ifstreamer is;
    
    is.open(fname);

    int maxsize;

    is.read(maxsize);
    for (int m = 0; m < maxsize; ++m) {
        ipoint2 key;
        max_vector_t& mapped = maxima[key];
        int vsize;

        is.read(key);
        is.read(vsize);
        mapped.reserve(vsize);
        for (int v = 0; v < vsize; ++v) {
            max_item_t item;

            is.read(item.pos);
            is.read(item.dist);
            mapped.push_back(item);
        }
    }
    is.read(regions);
    is.read(rcenters);

    is.close();
}
