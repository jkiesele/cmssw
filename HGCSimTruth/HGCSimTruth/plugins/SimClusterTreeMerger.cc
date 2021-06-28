#include <memory>
#include <vector>
#include <cstdlib>
#include <iostream>
#include <stack>
#include <unordered_map>
#include <map>
#include <sstream>
#include <utility>
#include <set>
#include <cmath>
#include <numeric>
using std::vector;
using std::map;
using std::unordered_map;
using std::pair;

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/Common/interface/Ref.h"

#include "SimDataFormats/Track/interface/SimTrack.h"
#include "SimDataFormats/Vertex/interface/SimVertex.h"
#include "SimDataFormats/Track/interface/SimTrackContainer.h"
#include "SimDataFormats/Vertex/interface/SimVertexContainer.h"
#include "SimDataFormats/CaloAnalysis/interface/SimClusterFwd.h"
#include "SimDataFormats/CaloAnalysis/interface/SimCluster.h"

#include "RecoLocalCalo/HGCalRecAlgos/interface/RecHitTools.h"

#include "DataFormats/Common/interface/Ptr.h"
#include "DataFormats/Common/interface/View.h"
#include "SimDataFormats/CaloHit/interface/PCaloHitContainer.h"
#include "SimDataFormats/CaloHit/interface/PCaloHit.h"
#include "DataFormats/DetId/interface/DetId.h"

#include "DataFormats/Common/interface/Association.h"
#include "DataFormats/Common/interface/AssociationMap.h"
#include "DataFormats/Common/interface/OneToManyWithQualityGeneric.h"

#include <iterator> // For std::forward_iterator_tag
#include <cstddef>  // For std::ptrdiff_t


//DIRTY PLOTTING - REMOVE LATER
#include "TCanvas.h"
#include "TH1F.h"
#include "TH2D.h"
#include "TEllipse.h"
#include <cmath>
//END DIRTY PLOTTING

#define EDM_ML_DEBUG
#define PI 3.14159265358979323846

#define MINOVERLAP 0.3
#define MINCIRCLERADIUS 1.

typedef edm::AssociationMap<edm::OneToManyWithQualityGeneric<
    SimClusterCollection, SimClusterCollection, float>> SimClusterToSimClusters;

/*
Yet another implementation of a 3D vector in CMSSW, but somehow
all the existing classes are missing needed functionality for the
merging algorithm.
*/
struct Vector3D{
    Vector3D() {}
    Vector3D(double x, double y, double z) : x_(x), y_(y), z_(z) {}
    Vector3D(const GlobalPoint& p) : x_(p.x()), y_(p.y()), z_(p.z()) {}
    Vector3D(const Vector3D& p) : x_(p.x_), y_(p.y_), z_(p.z_) {}
    ~Vector3D() {}
    Vector3D operator+(const Vector3D& other) const { return Vector3D(x_+other.x_, y_+other.y_, z_+other.z_); }
    Vector3D operator-(const Vector3D& other) const { return Vector3D(x_-other.x_, y_-other.y_, z_-other.z_); }
    Vector3D operator/(const double c){ return Vector3D(x_/c, y_/c, z_/c); }

    double norm(){ return std::sqrt(std::pow(x_,2) + std::pow(y_,2) + std::pow(z_,2)); }
    double dot(const Vector3D& other){ return x_ * other.x_ + y_ * other.y_ + z_ * other.z_; }

    double x_, y_, z_;
    friend std::ostream& operator<<(std::ostream& os, const Vector3D& v);
    };

std::ostream& operator<<(std::ostream& os, const Vector3D& v){
    os << "(" << v.x_ << ", " << v.y_ << ", " << v.z_ << ")";
    return os;
    }

Vector3D operator*(const double c, const Vector3D& p){ return Vector3D(c*p.x_, c*p.y_, c*p.z_); }
vector<Vector3D> operator+(const vector<Vector3D>& ps, const Vector3D& q){
    vector<Vector3D> out(ps.size());
    for (std::size_t i = 0; i < ps.size(); ++i){ 
        const Vector3D& p = ps[i];
        out[i] = p + q;
        }
    return out;
    }
vector<Vector3D> operator-(const vector<Vector3D>& ps, const Vector3D& q){
    vector<Vector3D> out(ps.size());
    for (std::size_t i = 0; i < ps.size(); ++i){ 
        const Vector3D& p = ps[i];
        out[i] = p - q;
        }
    return out;
    }

/* Needed for elementary 3D rotations */
struct RotMat3D{
    RotMat3D() {}
    ~RotMat3D() {}
    RotMat3D(
        double e11, double e12, double e13,
        double e21, double e22, double e23,
        double e31, double e32, double e33
        ) :
        e11_(e11), e12_(e12), e13_(e13),
        e21_(e21), e22_(e22), e23_(e23),
        e31_(e31), e32_(e32), e33_(e33)
        {}

    Vector3D dot(const Vector3D& p)const{
        return Vector3D(
            e11_*p.x_ + e12_*p.y_ + e13_*p.z_,
            e21_*p.x_ + e22_*p.y_ + e23_*p.z_,
            e31_*p.x_ + e32_*p.y_ + e33_*p.z_
            );
        }

    vector<Vector3D> dot(const vector<Vector3D>& ps)const{
        vector<Vector3D> out(ps.size());
        for (std::size_t i = 0; i < ps.size(); ++i){
            out[i] = dot(ps[i]);
            }
        return out;
        }

    RotMat3D dot(const RotMat3D& o)const{
        return RotMat3D(
            e11_*o.e11_ + e12_*o.e21_ + e13_*o.e31_, e11_*o.e12_ + e12_*o.e22_ + e13_*o.e32_, e11_*o.e13_ + e12_*o.e23_ + e13_*o.e33_,
            e21_*o.e11_ + e22_*o.e21_ + e23_*o.e31_, e21_*o.e12_ + e22_*o.e22_ + e23_*o.e32_, e21_*o.e13_ + e22_*o.e23_ + e23_*o.e33_,
            e31_*o.e11_ + e32_*o.e21_ + e33_*o.e31_, e31_*o.e12_ + e32_*o.e22_ + e33_*o.e32_, e31_*o.e13_ + e32_*o.e23_ + e33_*o.e33_
            );
        }

    RotMat3D transpose()const{
        return RotMat3D(
            e11_, e21_, e31_,
            e12_, e22_, e32_,
            e13_, e23_, e33_
            );
        }

    friend std::ostream& operator<<(std::ostream& os, const RotMat3D& m);

    double e11_, e12_, e13_,
          e21_, e22_, e23_,
          e31_, e32_, e33_;
    };

std::ostream& operator<<(std::ostream& os, const RotMat3D& m){
    os  << "[" << m.e11_ << ", " << m.e12_ << ", " << m.e13_ << "]\n"
        << "[" << m.e21_ << ", " << m.e22_ << ", " << m.e23_ << "]\n"
        << "[" << m.e31_ << ", " << m.e32_ << ", " << m.e33_ << "]";
    return os;
    }


struct Hit {
    Hit(double x, double y, double z, double t, double energy, int trackid, int layer, float radius, bool issilicon) :
        x_(x), y_(y), z_(z), t_(t), energy_(energy), trackid_(trackid),layer_(layer),radius_(radius),
        issilicon_(issilicon){}
    ~Hit() {}
    double x_;
    double y_;
    double z_;
    double t_;
    double energy_;
    int trackid_;
    int layer_;
    GlobalPoint gpoint(){return GlobalPoint(x_, y_, z_);}
    Vector3D vector3d(){return Vector3D(x_, y_, z_);}
    float radius_;
    bool issilicon_;
    };

/* Computes the 'average' position of a list of hits */
Vector3D hitcentroid(vector<Hit*>& hits){
    if (hits.size()==0) cms::Exception("SimMerging") << "Cannot compute hit centroid for 0 hits";
    else if (hits.size()==1) return Vector3D(hits[0]->x_, hits[0]->y_, hits[0]->z_);
    double summedEnergy = 0.;
    for(auto hit : hits) summedEnergy += hit->energy_;
    double center_x = 0., center_y = 0., center_z = 0.;
    for(auto hit : hits){
        double weight = hit->energy_/summedEnergy;
        center_x += weight * hit->x_;
        center_y += weight * hit->y_;
        center_z += weight * hit->z_;
        }
    return Vector3D(center_x, center_y, center_z);
    }

/* Cumulative sum of a vector */
template<typename T>
vector<T> cumsum(const vector<T>& input){
    T s = 0.;
    vector<T> cumulative_sum(input.size());
    for (std::size_t i = 0; i < input.size(); ++i) {
        s += input[i];
        cumulative_sum[i] = s;
        }
    return cumulative_sum;
    }

/*
Returns a permutation of a sorted vector, which is applicable on another vector
See https://stackoverflow.com/a/17074810
*/
template <typename T, typename Compare>
vector<std::size_t> argsort(
    const vector<T>& vec,
    Compare& compare)
    {
        vector<std::size_t> p(vec.size());
        std::iota(p.begin(), p.end(), 0);
        std::sort(p.begin(), p.end(),
            [&](std::size_t i, std::size_t j){ return compare(vec[i], vec[j]); });
        return p;
        }
/* Same but using a default "<" comparison */
template <typename T>
vector<std::size_t> argsort(const vector<T>& vec){
    vector<std::size_t> p(vec.size());
    std::iota(p.begin(), p.end(), 0);
    std::sort(
        p.begin(), p.end(),
        [&](std::size_t i, std::size_t j){ return vec[i] < vec[j]; }
        );
    return p;
    }
/* Function to apply the argsort returned from argsort(), returning a copy */
template <typename T>
std::vector<T> apply_argsort(
    const std::vector<T>& vec,
    const std::vector<std::size_t>& p)
{
    std::vector<T> sorted_vec(vec.size());
    std::transform(p.begin(), p.end(), sorted_vec.begin(),
        [&](std::size_t i){ return vec[i]; });
    return sorted_vec;
}
/* Function to apply the argsort returned from argsort() in place */
template <typename T>
void apply_argsort_in_place(
    std::vector<T>& vec,
    const std::vector<std::size_t>& p)
{
    std::vector<bool> done(vec.size());
    for (std::size_t i = 0; i < vec.size(); ++i)
    {
        if (done[i])
        {
            continue;
        }
        done[i] = true;
        std::size_t prev_j = i;
        std::size_t j = p[i];
        while (i != j)
        {
            std::swap(vec[prev_j], vec[j]);
            done[j] = true;
            prev_j = j;
            j = p[j];
        }
    }
}

class Node {
    public:
        Node() : 
            trackid_(0), pdgid_(0), initial_energy_(0.), parent_(nullptr)
            {}
        Node(const SimTrack& track,vector<double> merging_threshold_transv,vector<double> merging_threshold_longitud) {
            trackid_ = track.trackId();
            initial_energy_ = track.momentum().E();
            pdgid_ = track.type();
            is_hadron_ = pdgid_ > 111;
            final_z_ = track.trackerSurfacePosition().z();
            merged_trackids_.push_back(trackid_);
            crossed_boundary_ = track.crossedBoundary();
            if (crossed_boundary_){
                boundary_momentum_ = track.getMomentumAtBoundary();
                boundary_position_ = Vector3D(
                    track.getPositionAtBoundary().x(),
                    track.getPositionAtBoundary().y(),
                    track.getPositionAtBoundary().z()
                    );
                boundary_time_ = track.getPositionAtBoundary().t();
                centroid_=boundary_position_;
                }
            parent_=0;
            circle_radius_=0.1;
            hadhits_=false;
            isdense_=false;
            }
        ~Node() {}

        /* Number of quantities that depend on the hits */
        void calculate_shower_variables(){
            if (nhits() == 0) return;

            if(!isLeaf()){
               // circle_radius_=0;
               // std::cout << "calculating shower variables for non leaf node?" << std::endl;//DEBUG
               // return;
            }
            //DIRTY HACKS HERE
            auto hitcp = hits_;
            //find lower layer number
            int minlayer=100;
            double highestenergy=0;
            for(const auto h:hits_){
                if(h->energy_>highestenergy)
                    highestenergy=h->energy_;
            }
            double energycut = std::min(0.00005, highestenergy);
            for(const auto h:hits_){
                if(h->energy_ < energycut)
                    continue;
                if(h->layer_<minlayer)
                    minlayer = h->layer_;
            }
            //now select only the first 5 layers
            vector<Hit*> newhits;
            for(const auto h:hits_){
                if(h->layer_<minlayer+1){
                    newhits.push_back(h);
                }
            }
            hits_=newhits;

            if(hits_.size()<7)
                return;

            centroid_ = ( ::hitcentroid(hits_));

            //is it the axis? this only works for straight lines..
            //pca like would be useful here FIXME
            //but maybe taking this as first approx, and then calculate correction rather than full PCA?
            auto axis_  = ((centroid_-boundary_position_) / (centroid_-boundary_position_).norm());


            // Compute the transverse distances from hits to the shower axis
            vector<double> d_to_axis;
            vector<double> d_along_axis;
            vector<double> energies;
            double total_energy = 0.;
            for(auto hit : hits_){
                Vector3D hit_pos = hit->vector3d() - boundary_position_;
                Vector3D projection_along_axis = hit_pos.dot(axis_) * axis_;
                d_along_axis.push_back(projection_along_axis.norm());
                //NEW!!
                //DIRTY HACK FIXME just for testing DEBUG
                double thisd= 2.*hit->radius_ + (hit_pos - projection_along_axis).norm();
                d_to_axis.push_back(thisd);
                //notice 2* radius gets added, and is therefore kind of a lower limit
                energies.push_back(hit->energy_);
                total_energy += hit->energy_;
            }

            vector<std::size_t> order = argsort(d_to_axis);
            apply_argsort_in_place(d_to_axis, order);
            vector<double> cumsum_energies_to_axis = cumsum(apply_argsort(energies, order));
            auto ordered_hits = apply_argsort(hits_, order);

            // Find the energy containment radii
            // vector<double> thresholds = { .3, .75, .85 };

            //FIXME, we need to actually use the right radii here!
            //but shouldn't change the execution time a lot

            isdense_=false;
            for (std::size_t i = 0; i < cumsum_energies_to_axis.size()-1; ++i) {
                double cumsum_this = cumsum_energies_to_axis[i] / total_energy;
                double cumsum_next = cumsum_energies_to_axis[i+1] / total_energy;
                double threshold = 0.1;//scale later
                circle_radius_ = d_to_axis[i+1];//DEBUG
                if(d_to_axis[i] < 5.*ordered_hits[i]->radius_){//is there any sort of dense core?
                    isdense_=true;
                }
                if (cumsum_next > threshold)
                    break;

            }
            if(circle_radius_<MINCIRCLERADIUS)
                circle_radius_=MINCIRCLERADIUS;//


        }

        /* Standard depth-first-search tree traversal as an iterator */
        struct Iterator {
            using iterator_category = std::forward_iterator_tag;
            using difference_type   = std::ptrdiff_t;
            using value_type        = Node;
            using pointer           = Node*;  // or also value_type*
            using reference         = Node&;  // or also value_type&

            Iterator(pointer ptr, bool verbose=false) :
                m_ptr(ptr), root(ptr), depth_(0), verbose_(verbose) {}

            reference operator*() const { return *m_ptr; }
            pointer operator->() { return m_ptr; }

            Iterator& operator++() {
                if (m_ptr->hasChildren()){
                    if (verbose_) edm::LogVerbatim("SimMerging")
                        << "Track " << m_ptr->trackid_
                        << ": Going to first child " << m_ptr->children_[0]->trackid_
                        ;
                    continuation_.push(m_ptr);
                    m_ptr = m_ptr->children_[0];
                    depth_++;
                    }
                else {
                    if (verbose_) edm::LogVerbatim("SimMerging")
                        << "Track " << m_ptr->trackid_
                        << ": No children, going to next sibling"
                        ;
                    while(true){
                        if (m_ptr == root){
                            if (verbose_) edm::LogVerbatim("SimMerging") << "Back at the root of the iterator; quiting";
                            m_ptr = nullptr;
                            break;
                            }
                        else if (m_ptr->hasNextSibling()){
                            m_ptr = m_ptr->nextSibling();
                            if (verbose_) edm::LogVerbatim("SimMerging") << "Has sibling; going to " << m_ptr->trackid_;
                            break;
                            }
                        if (verbose_) edm::LogVerbatim("SimMerging") << "Has no sibling; proceed popping stack";
                        m_ptr = continuation_.top();
                        continuation_.pop();
                        depth_--;
                        if (verbose_) edm::LogVerbatim("SimMerging") << "Popped back to track " << m_ptr->trackid_;
                        }
                    }
                return *this;
                }
            // Postfix increment
            Iterator operator++(int) { Iterator tmp = *this; ++(*this); return tmp; }
            const int depth(){return depth_;}
            friend bool operator== (const Iterator& a, const Iterator& b) { return a.m_ptr == b.m_ptr; };
            friend bool operator!= (const Iterator& a, const Iterator& b) { return a.m_ptr != b.m_ptr; };  
            private:
                pointer m_ptr;
                pointer root;
                int depth_;
                bool verbose_;
                std::stack<pointer> continuation_;
            };
        Iterator begin(bool verbose=false) { return Iterator(this, verbose); }
        Iterator end() { return Iterator(nullptr); }

        /* Traverses upwards */
        struct IteratorUp {
            using iterator_category = std::forward_iterator_tag;
            using difference_type   = std::ptrdiff_t;
            using value_type        = Node;
            using pointer           = Node*;  // or also value_type*
            using reference         = Node&;  // or also value_type&
            IteratorUp(pointer ptr) : m_ptr(ptr) {}
            reference operator*() const { return *m_ptr; }
            pointer operator->() { return m_ptr; }
            IteratorUp& operator++() {
                m_ptr = (m_ptr->hasParent()) ? m_ptr->parent_ : nullptr ;
                return *this;
                }
            IteratorUp operator++(int) { IteratorUp tmp = *this; ++(*this); return tmp; }
            friend bool operator== (const IteratorUp& a, const IteratorUp& b) { return a.m_ptr == b.m_ptr; };
            friend bool operator!= (const IteratorUp& a, const IteratorUp& b) { return a.m_ptr != b.m_ptr; };  
            private: pointer m_ptr;
            };
        IteratorUp begin_up() { return IteratorUp(this); }
        IteratorUp end_up() { return IteratorUp(nullptr); }

        void setParent(Node* parent) {parent_ = parent;}
        void addChild(Node* child) {children_.push_back(child);}
        void addHit(Hit* hit) {hits_.push_back(hit);hadhits_=true;}
        int nhits()const{ return hits_.size(); }
        bool hadHits()const{ return hadhits_; }
        void clearHits(){hits_.clear();}
        bool isLeaf(){ return children_.empty(); }
        bool hasChildren(){ return !(children_.empty()); }
        bool hasParent(){ return parent_ != nullptr; }

        bool hasNextSibling(){
            if (!parent_){
                // There is no parent
                return false;
                }
            else if (parent_->children_.back() == this){
                // This was the last child
                return false;
                }
            // else if ( (int)(parent_->children_.size()-1) == iChild_){
            //     // This was the last child
            //     return false;
            //     }
            return true;
            }

        Node* nextSibling(){
            if (hasNextSibling()){
                vector<Node*>::iterator sibling = std::find(
                    parent_->children_.begin(), parent_->children_.end(), this
                    );
                sibling++; // advance once
                if (sibling == parent_->children_.end())
                    // This shouldn't happen
                    return nullptr;
                return *sibling;
                }
            return nullptr;
            }

        /* Traverses tree and builds string representation */
        std::string stringrep(){
            std::stringstream ss;
            int nTracks = 0;
            int nHits = 0;
            for (Node::Iterator it = begin(); it != end(); it++){
                Node& node = (*it);
                for (int i = 0; i < it.depth(); ++i){
                    ss << "--";
                    }
                ss
                    << "Track " << node.trackid_
                    << " (" << node.nhits() << " hits)"
                    << "\n";
                nTracks++;
                nHits += node.nhits();
                }
            ss << "In total " << nTracks << " tracks with " << nHits << " hits";
            return ss.str();
            }

        /* A node is a 'leaf parent' if it has children, and all those children are leafs  */
        bool isLeafParent(){
            // A leaf itself is not a leaf parent
            if (isLeaf()) return false;
            for(auto child : children_){
                if (child->hasChildren()) return false;
                }
            return true;
            }

        /*
        Flips all basic z quantities by -1
        Does not flip the shower properties; calculate_shower_variables must be rerun
        */
        void flipz(){
            final_z_ *= -1.;
            if(crossed_boundary_)
                boundary_position_.z_ *= -1.;
            for(auto & h : hits_)
                h->z_ *= -1;
        //?    boundary_momentum_.SetPz(-boundary_momentum_.z());
        }


        void absorbNode(Node *rhs){

            for (auto trackid : rhs->merged_trackids_){
                merged_trackids_.push_back(trackid);
            }

            for(auto child : rhs->children_){
                addChild(child);
                child->setParent(this);
            }
            rhs->children_.clear();

            //non dense nodes can *be* merged but cannot contribute to area after first merge anymore
            if(!isdense_){
                circle_radius_=0.1;//set low
            }
            if(!rhs->isdense_){
                rhs->circle_radius_=0.1;
            }

            hadhits_ |= rhs->hadhits_;
        }

        bool crossed_boundary_, is_hadron_;
        int trackid_, pdgid_;
        double initial_energy_, final_z_;

        //use children pointers for all circles etc.
        //use hits only temporarily, don't store them (second iteration)
        //that makes:
        math::XYZTLorentzVectorF boundary_momentum_;
        Vector3D boundary_position_;
        double boundary_time_;//this could contribute include time
        double circle_radius_; //just calculate once, using PDGID

        Node * parent_;
        vector<Node*> children_;
        vector<int> merged_trackids_;
        bool hadhits_;
        vector<Hit*> hits_;
        bool isdense_;

        //just for debugging/dev
        Vector3D centroid_;

    };

/* Finds a track by trackid in a tree */
Node* find_in_tree(Node* root, int trackid){
    for(auto& track : *root){
        if (track.trackid_ == trackid){
            return &track;
            }
        }
    throw cms::Exception("Unknown")
        << "Track id " << trackid << " is not in the tree";
    }


double calculate_circle_overlap(const double& d, const double& r1,
        const double& r2){
//from https://www.xarg.org/2016/07/calculate-the-intersection-area-of-two-circles/
    if (d < r1 + r2) {

        double a = r1 * r1;
        double b = r2 * r2;

        double x = (a - b + d * d) / (2 * d);
        double z = x * x;
        double y = sqrt(a - z);

        if (d <= std::abs(r2 - r1)) {
            return PI * std::min(a, b);
        }
        return a * asin(y / r1) + b * asin(y / r2) - y * (x + sqrt(z + b - a));
    }
    return 0;
}
double calculate_shower_overlap(Node* t1, Node* t2){
    //now overlap with all children here; just one-by-one, not of union
    //which will be more complicated

    //maybe making this an explicit 'root' check plus childred check is better than copying vectors around
    vector<Node*> nodes_a = t1->children_;
    nodes_a.push_back(t1);
    vector<Node*> nodes_b = t2->children_;
    nodes_b.push_back(t2);

    double max_overlap=0;

    for(const auto & na: nodes_a){
        if(!na->circle_radius_)
            continue;
        for(const auto& nb: nodes_b){
            if(!nb->circle_radius_)//non leafs
                continue;

            //we could reintroduce the cachedOverlap thing here

            //quick check if there can be any overlap
            auto distvec = na->boundary_position_ - nb->boundary_position_;
            double distsq = distvec.dot(distvec); //x^2 OPT possible (just take x and y)
            double rsum = na->circle_radius_+nb->circle_radius_;

            if(rsum*rsum<distsq)
                continue; //can't overlap


            //JUST SIMPLE. either hits the others centre: this
            //works because the sensor radius is part of the circles
            auto maxrad = std::max(na->circle_radius_,nb->circle_radius_);
            if(distsq<maxrad*maxrad)
                max_overlap=1;
            continue;

            //maybe this doesn't even need to be area but just a distance measure!
            // (OPT possible) because overlap area is monotone function of that

            double overlap = calculate_circle_overlap(std::sqrt(distsq),
                    na->circle_radius_,nb->circle_radius_);
            //normalise overlap here
            double minrad = std::min(na->circle_radius_, nb->circle_radius_);

            overlap /= PI * minrad*minrad;

            //std::cout << "rel overlap " << overlap << std::endl; //DEBUG
            //do something with boundary_time_ here if you want

            if(overlap>max_overlap)
                max_overlap=overlap;

        }
    }
    return max_overlap;
}



/* Remove a node from its parent's children vector */
void break_from_parent(Node* node){
    if (!(node->hasParent())) cms::Exception("SimMerging") << "Cannot remove root node";
    vector<Node*>& parents_children = node->parent_->children_;
    // erase-remove idiom: https://en.wikipedia.org/wiki/Erase%E2%80%93remove_idiom#Example
    parents_children.erase(
        std::remove(parents_children.begin(), parents_children.end(), node),
        parents_children.end()
        );
    }

/* Breaks node from parent but moves its children to the children of the parent */
void remove_intermediate_node(Node* node){
    Node* parent = node->parent_;
    break_from_parent(node);
    // Move children of the now-removed node to its parent
    for (auto child : node->children_){
        parent->addChild(child);
        child->setParent(parent);
        }
    }

// _______________________________________________
// Some functions for traversal by recursion
// These first build the whole traversal in a vector
// (as pointers, so memory usage is not too bad)

/* Does depth first search traversal by using recursion, but not as an iterator */
void _dfs_recursion(
        Node* node,
        vector<pair<Node*, int>>& returnable,
        int depth
        )
    {
    returnable.push_back(std::make_pair(node, depth));
    for (auto child : node->children_){
        _dfs_recursion(child, returnable, depth+1);
        }
    }

/*
Wrapper around the recursive version that only takes a node as input.
Returns a vector of pair<node, depth>.
Useful if you want to keep the whole traversal in memory; usually you
will want to use the iterator version of the Node class.
*/
vector<pair<Node*, int>> dfs(Node* root){
    vector<pair<Node*, int>> returnable;
    _dfs_recursion(root, returnable, 0);
    return returnable;
    }

/* String representation of dfs traversal (keeps whole traversal in memory) */
std::string dfs_stringrep(Node* root){
    std::stringstream ss;
    for (auto node_depth_pair : dfs(root)){
        for (int i = 0; i < node_depth_pair.second; ++i){
            ss << "--";
            }
        ss
            << "Track " << node_depth_pair.first->trackid_
            << " (" << node_depth_pair.first->nhits() << " hits)"
            << "\n";
        }
    return ss.str();
    }

/* Remove single-child no-hit tracks (i.e. intermediate tracks) */
void trim_tree(Node* root){
    // Traverse once and note all tracks that should be kept:
    // Either a track that has hits, or an ancestor thereof
    std::set<int> trackids_with_hits_or_parents_thereof;
    for (auto& node : *root){
        if (!(node.hadHits())) continue;
        // Iterate upwards and save in the set
        for (auto it=node.begin_up(); it!=node.end_up(); it++){
            trackids_with_hits_or_parents_thereof.insert(it->trackid_);
            }
        }
    // Now remove all nodes not in the set
    // We'll be modifying parent-child relationships mid-loop, so we have to be a little
    // careful
    auto it=root->begin();
    while(it!=root->end()){
        Node& node = (*it);
        if (!(trackids_with_hits_or_parents_thereof.count(node.trackid_))){
            // First remove children so the iterator will go to the next sibling
            node.children_.clear();
            // Advance to next sibling (or further up the chain)
            it++;
            // Then break from parent (if doing this before advancing the order gets messed up)
            break_from_parent(&node);
            }
        else{
            it++;
            }
        }
    // Second trimming step: Remove 'intermediate' tracks
    // (i.e. tracks with no hits, 1 child, and 1 parent)
    // In this case it's easier to put the whole traversal in memory first,
    // and avoid modifying relationships mid-loop
    for (auto node_depth_pair : dfs(root)){
        Node* node = node_depth_pair.first;
        if (node->hasParent() && (node->children_.size()==1) && !(node->hadHits())){
            remove_intermediate_node(node);
            }
        }
    }


/* Compute a distance measure between two nodes: now simply distance between the hit centroids */
//not needed

bool merge_leafparent(Node* leafparent, double min_overlap){
    edm::LogVerbatim("SimMerging") << "  Merging leafparent " << leafparent->trackid_;
    bool didUpdate = false;
    // Copy list of potentially mergeable nodes
    vector<Node*> mergeable = leafparent->children_;
    leafparent->children_.clear();
    // Parent itself can be mergeable, if it has hits and is not a root
    if (leafparent->hasParent() && leafparent->hadHits()) mergeable.push_back(leafparent);
    // Start merging
    while(true){
        bool break_early = false;
        bool didUpdateThisIteration = false;
        double current_max_overlap = min_overlap;
        pair<Node*,Node*> pairToMerge;
        // Compute all distances between clusters
        int nMergeable = mergeable.size();
        for (int i = 0; i < nMergeable; ++i){
            Node* left = mergeable[i];
            for (int j = i+1; j < nMergeable; ++j){
                Node* right = mergeable[j];
                double overlap = calculate_shower_overlap(left,right);

                if (overlap > current_max_overlap){
                    current_max_overlap = overlap;
                    pairToMerge = (left->initial_energy_ > right->initial_energy_) ?
                            std::make_pair(left, right) : std::make_pair(right, left);
                    didUpdate = true;
                    didUpdateThisIteration = true;
                    if (overlap == 1.){
                        // It can't get higher than 1 anyway
                        break_early = true;
                        break;
                    }
                }

            }
            if (break_early) break;
        }
        if (!didUpdateThisIteration) break; // Nothing to merge this iteration
        // Now do the merging

        ///

        edm::LogVerbatim("SimMerging")
            << "    Merging " << pairToMerge.second->trackid_
            << " into " << pairToMerge.first->trackid_
            ;
        // Bookkeep that the track (and any previously merged tracks) is merged in

        pairToMerge.first->absorbNode(pairToMerge.second);
        break_from_parent(pairToMerge.second);
        mergeable.erase(
            std::remove(mergeable.begin(), mergeable.end(), pairToMerge.second),
            mergeable.end()
            );
        // Do NOT recompute the shower variables for newly merged node, now that it has more hits

        //add to first node
        //pairToMerge.first->calculate_shower_variables();


        ///

    }
    // Make a string representation of the mergeable nodes for debugging
    std::string mergeableStr = "";
    if(mergeable.size()){
        std::stringstream ss;
        for(auto child : mergeable) ss << child->trackid_ << ", ";
        mergeableStr = ss.str();
        mergeableStr.pop_back(); mergeableStr.pop_back(); // Remove trailing comma
        }
    // All possible merging now done;
    // Next steps depend on whether the passed node was a root
    if (!(leafparent->hasParent())){
        leafparent->children_ = mergeable;
        if(didUpdate) {
            // Simply overwrite with the merged nodes
            edm::LogVerbatim("SimMerging")
                << "    Root " << leafparent->trackid_
                << " is set to have the following children: "
                << mergeableStr;
            }
        else{
            edm::LogVerbatim("SimMerging")
                << "    Root " << leafparent->trackid_
                << ": no further merging possible";
            }
        return didUpdate;
        }
    else {
        // Special case: If the leafparent had no hits (and was thus not included as
        // a mergeable node), AND all nodes were merged into one cluster, assign the 
        // pdgid of the leafparent to the remaining node
        if(
            !(leafparent->hadHits())
            && mergeable.size()==1
            && mergeable[0]->pdgid_!=leafparent->pdgid_
            ){
            edm::LogVerbatim("SimMerging")
                << "    Using leafparent pdgid " << leafparent->pdgid_
                << " for track " << mergeable[0]->trackid_
                << " (rather than " << mergeable[0]->pdgid_
                << ") since all nodes were merged into one";
            mergeable[0]->pdgid_ = leafparent->pdgid_;
            }
        // Replace the node in the parent's children list with all merged nodes
        Node* parent = leafparent->parent_;
        edm::LogVerbatim("SimMerging")
            << "    Adding the following children to parent " << parent->trackid_
            << ": " << mergeableStr;
        break_from_parent(leafparent);
        for(auto child : mergeable){
            parent->addChild(child);
            child->setParent(parent);
            }
        return true;
        }
    }

void merging_algo(Node* root){
    int iIteration = -1;
    bool didUpdate = true;
    while(didUpdate){
        iIteration++;
        edm::LogVerbatim("SimMerging") << "Iteration " << iIteration;
        // Build list of leaf parents in memory
        vector<Node*> leafparents;
        for (auto& node : *root){
            if (!(node.isLeafParent())) continue;
            leafparents.push_back(&node);
            }
        for (auto node : leafparents){
            didUpdate = merge_leafparent(node, MINOVERLAP);
            }
        }
    edm::LogVerbatim("SimMerging") << "Done after iteration " << iIteration;
    }

// _______________________________________________


class simmerger : public edm::stream::EDProducer<> {
    public:
        explicit simmerger(const edm::ParameterSet&);
        ~simmerger() {}
        SimCluster mergedSimClusterFromTrackIds(std::vector<int>& trackIds, 
            const edm::Association<SimClusterCollection>& simTrackToSimCluster);
    private:
        virtual void produce(edm::Event&, const edm::EventSetup&) override;
        void beginRun(const edm::Run&, const edm::EventSetup&) override {}
        hgcal::RecHitTools hgcalRecHitToolInstance_ ;
        edm::EDGetTokenT<edm::View<PCaloHit>> hgcalEEHitsToken_;
        edm::EDGetTokenT<edm::View<PCaloHit>> hgcalHEfrontHitsToken_;
        edm::EDGetTokenT<edm::View<PCaloHit>> hgcalHEbackHitsToken_;
        edm::EDGetTokenT<edm::SimTrackContainer> tokenSimTracks;
        edm::EDGetTokenT<edm::SimVertexContainer> tokenSimVertices;
        edm::EDGetTokenT<SimClusterCollection> simClustersToken_;
        edm::EDGetTokenT<edm::Association<SimClusterCollection>> simTrackToSimClusterToken_;
        unordered_map<unsigned int, SimTrackRef> trackIdToTrackRef_;
        vector<double> mergeThresholdsTransv_,mergeThresholdsLongitud_;
    };


simmerger::simmerger(const edm::ParameterSet& iConfig) :
    hgcalEEHitsToken_(consumes<edm::View<PCaloHit>>(edm::InputTag("g4SimHits", "HGCHitsEE"))),
    hgcalHEfrontHitsToken_(consumes<edm::View<PCaloHit>>(edm::InputTag("g4SimHits", "HGCHitsHEfront"))),
    hgcalHEbackHitsToken_(consumes<edm::View<PCaloHit>>(edm::InputTag("g4SimHits", "HGCHitsHEback"))),
    tokenSimTracks(consumes<edm::SimTrackContainer>(edm::InputTag("g4SimHits"))),
    tokenSimVertices(consumes<edm::SimVertexContainer>(edm::InputTag("g4SimHits"))),
    simClustersToken_(consumes<SimClusterCollection>(edm::InputTag("mix:MergedCaloTruth"))),
    simTrackToSimClusterToken_(consumes<edm::Association<SimClusterCollection>>(edm::InputTag("mix:simTrackToSimCluster"))),
    mergeThresholdsTransv_(iConfig.getParameter<vector<double > > ( "MergeTheresholdsTransv" )),
    mergeThresholdsLongitud_(iConfig.getParameter<vector<double > > ( "MergeTheresholdsLongitud" ))
    {
    produces<SimClusterCollection>();
    produces<edm::Association<SimClusterCollection>>();
    produces<SimClusterToSimClusters>();
    }

SimCluster simmerger::mergedSimClusterFromTrackIds(std::vector<int>& trackIds, 
    const edm::Association<SimClusterCollection>& simTrackToSimCluster) {
    SimCluster sc; 
    for (auto tid : trackIds) {
        if (trackIdToTrackRef_.find(tid) == trackIdToTrackRef_.end())
            throw cms::Exception("SimClusterTreeMerger") << "Failed to find a trackId in the TrackMap.";
        sc += *simTrackToSimCluster[trackIdToTrackRef_[tid]];
    }
    return sc;
}

static Vector3D project_to_boundary(const Vector3D& toprop, const Vector3D& bpos){

    SimCluster::LorentzVector vec (toprop.x_,toprop.y_, toprop.z_, 0.);
    auto r = bpos.z_*tan(vec.theta());

    SimCluster::PolarLorentzVector pvec(r, vec.Eta(), vec.Phi(), 0);
    return Vector3D(pvec.x(),pvec.y(),pvec.z());
}

void simmerger::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {  
    edm::ESHandle<CaloGeometry> geom;
    iSetup.get<CaloGeometryRecord>().get(geom);
    hgcalRecHitToolInstance_.setGeometry(*geom);
    trackIdToTrackRef_.clear();

    auto output = std::make_unique<SimClusterCollection>();

    // Create Hit instances
    vector<Hit> hits;
    vector<edm::EDGetTokenT<edm::View<PCaloHit>>> tokens = {
        hgcalEEHitsToken_,
        hgcalHEfrontHitsToken_,
        hgcalHEbackHitsToken_
        };
    std::set<int> trackids_with_hits;
    for (edm::EDGetTokenT<edm::View<PCaloHit>> token : tokens ) {
        edm::Handle< edm::View<PCaloHit> > handle;
        iEvent.getByToken(token, handle);
        for (auto const & hit : handle->ptrs() ) {
            DetId id = hit->id();
            GlobalPoint position = hgcalRecHitToolInstance_.getPosition(id);
            float radius = 2;//DEBUG FIXME needs scintillator radii
            if(hgcalRecHitToolInstance_.isSilicon(id))
                radius=hgcalRecHitToolInstance_.getRadiusToSide(id);

            if(radius<0.5)
                radius=0.5;
            hits.push_back(Hit(
                position.x(), position.y(), position.z(),
                hit->time(), hit->energy(), hit->geantTrackId(),
                hgcalRecHitToolInstance_.getLayer(id),
                radius,
                hgcalRecHitToolInstance_.isSilicon(id)
                ));
            trackids_with_hits.insert(hit->geantTrackId());
            }
        }

    // Build the tree
    edm::Handle<edm::SimTrackContainer> handleSimTracks;
    iEvent.getByLabel("g4SimHits", handleSimTracks);
    edm::Handle<edm::SimVertexContainer> handleSimVertices;
    iEvent.getByLabel("g4SimHits", handleSimVertices);

    edm::LogVerbatim("SimMerging") << "Building map";
    unordered_map<int, Node> trackid_to_node;
    for(size_t i = 0; i < handleSimTracks->size(); i++){
        SimTrackRef track(handleSimTracks, i);
        trackid_to_node.emplace(track->trackId(), Node(*track, mergeThresholdsTransv_,mergeThresholdsLongitud_));
        trackIdToTrackRef_[track->trackId()] = track;
	}

    edm::LogVerbatim("SimMerging") << "Adding hits to nodes";
    for (auto& hit : hits){
        trackid_to_node[hit.trackid_].addHit(&hit);
        }

    edm::LogVerbatim("SimMerging") << "Building tree";
    Node* root(new Node());
    for(const auto& track : *handleSimTracks){
        int trackid = track.trackId();
        Node* node = &(trackid_to_node[trackid]);
        // Have to get parent info via the SimVertex
        SimVertex vertex = handleSimVertices.product()->at(track.vertIndex());
        bool hasParent = !(vertex.noParent());
        // Build the tree
        if (hasParent){
            int parentid = vertex.parentIndex();
            edm::LogVerbatim("SimMerging")
                << "Setting parent->child relationship: "
                << parentid << " -> " << trackid
                ;
            auto it = trackid_to_node.find(parentid);
            if (it != trackid_to_node.end()){
                Node* parent = &(it->second);
                node->setParent(parent);
                parent->addChild(node);
                }
            else{
                throw cms::Exception("Unknown")
                    << "Track id " << parentid
                    << " is not in the map"
                    ;
                }
            }
        else{
            edm::LogVerbatim("SimMerging") << "Found parentless particle: " << node->trackid_;
            root->addChild(node);
            node->setParent(root);
            }
        }

#ifdef EDM_ML_DEBUG
    edm::LogVerbatim("SimMerging") << "Printing root " << root->trackid_;
    edm::LogVerbatim("SimMerging") << root->stringrep() << "\n";
    edm::LogVerbatim("SimMerging") << "Trimming tree...";
#endif

    trim_tree(root);

#ifdef EDM_ML_DEBUG
    edm::LogVerbatim("SimMerging") << "Printing root " << root->trackid_ << " after trimming";
    edm::LogVerbatim("SimMerging") << root->stringrep() << "\n";
#endif

#ifdef EDM_ML_DEBUG
    edm::LogVerbatim("SimMerging") << "Splitting into positive and negative endcap";
#endif


    //Jan comment: presumably that's where we can split by particle
    Node* pos(new Node());
    Node* neg(new Node());
    for(auto child: root->children_){
        if(child->final_z_ < 0.){
            child->setParent(neg);
            neg->addChild(child);
            }
        else{
            child->setParent(pos);
            pos->addChild(child);
            }
        }

#ifdef EDM_ML_DEBUG
    edm::LogVerbatim("SimMerging") << "Flipping all z's in negative endcap for merging algorithm";
#endif



#ifdef EDM_ML_DEBUG
    edm::LogVerbatim("SimMerging") << "Printing neg " << neg->trackid_ << " after merging_algo";
    edm::LogVerbatim("SimMerging") << neg->stringrep() << "\n";
#endif


    //DIRTY plot

    std::cout << "pos side " << std::endl; //DEBUG


    static int counter=0;
    auto s= std::to_string(counter);
    counter++;


    TCanvas cv;
    TH2D h("","",100,-200,200,100,-200,200);
    std::vector<TEllipse*> allell;

    for (auto& node : *pos){
        node.calculate_shower_variables();

        //plotting here, can be removed later
        if(node.hadHits()){//plot one bigone
            auto bpos = node.boundary_position_;
            for(const auto& hit:node.hits_){
                //simple eta phi projection
                auto projvec = project_to_boundary(Vector3D(hit->x_,hit->y_, hit->z_),bpos);
                h.Fill(projvec.x_,projvec.y_);

            }
            auto projvec = project_to_boundary(node.centroid_,node.boundary_position_);

            auto ell = new TEllipse(projvec.x_,
                    projvec.y_, node.circle_radius_, node.circle_radius_);
            if(!node.isdense_)
                ell->SetLineColor(kRed);
            ell->SetFillStyle(0);
            allell.push_back(ell);
        }


        node.clearHits();

    }
    h.Draw("colz");
    cv.SetLogz();
    for(auto& el:allell)
        el->Draw("same");
    cv.Print((TString)"test_all"+s+".pdf");
    for(auto e:allell)
        delete e;



    ///DIRTY PLOT END

    merging_algo(pos);

    std::cout << "neg side " << std::endl; //DEBUG
    for (auto& node : *neg){
        node.flipz();
        node.calculate_shower_variables();
        node.clearHits();
    }
    merging_algo(neg);

    // Fill the output; the clusters are the remaining nodes (except the root)
    edm::Handle<edm::Association<SimClusterCollection>> simTrackToSimClusterHandle;
    iEvent.getByToken(simTrackToSimClusterToken_, simTrackToSimClusterHandle);

    edm::Handle<SimClusterCollection> simClusterHandle;
    iEvent.getByToken(simClustersToken_, simClusterHandle);

	auto allChildren = pos->children_;
	allChildren.insert(allChildren.begin(), neg->children_.begin(), neg->children_.end());

    size_t i = 0;
    std::vector<int> mergedIndices(simClusterHandle->size(), 0);
    std::vector<std::vector<SimClusterRef>> groups;
    for(auto cluster : allChildren) {
        SimCluster sc; 
        std::vector<SimClusterRef> group;
        for (auto tid : cluster->merged_trackids_) {
            if (trackIdToTrackRef_.find(tid) == trackIdToTrackRef_.end())
                throw cms::Exception("SimClusterTreeMerger") << "Failed to find a trackId in the TrackMap.";
            const auto& unmerged = (*simTrackToSimClusterHandle)[trackIdToTrackRef_[tid]];
            mergedIndices.at(unmerged.key()) = i;
            sc += *unmerged;
            group.push_back(unmerged);
        }
        std::sort(std::begin(group), std::end(group), 
            [](SimClusterRef a, SimClusterRef b) { return a->impactMomentum().energy() > b->impactMomentum().energy(); });
        i++;
		sc.setPdgId(cluster->pdgid_);
        output->push_back(sc);
        groups.push_back(group);
    }

    const auto& mergedSCHandle = iEvent.put(std::move(output));

    auto assocMap = std::make_unique<SimClusterToSimClusters>(mergedSCHandle, simClusterHandle);

    for (size_t i = 0; i < groups.size(); i++) {
        SimClusterRef msc(mergedSCHandle, i);
        float mergedE = msc->impactMomentum().energy();
        for (auto sc : groups.at(i)) {
            float scE = sc->impactMomentum().energy();
            assocMap->insert(msc, std::make_pair(sc, mergedE/scE));
        }
    }
    iEvent.put(std::move(assocMap));

    auto assoc = std::make_unique<edm::Association<SimClusterCollection>>(mergedSCHandle);
    edm::Association<SimClusterCollection>::Filler filler(*assoc);
    filler.insert(simClusterHandle, mergedIndices.begin(), mergedIndices.end());
    filler.fill();
    iEvent.put(std::move(assoc));

    delete root;
    delete pos;
    delete neg;
    }

DEFINE_FWK_MODULE(simmerger);
