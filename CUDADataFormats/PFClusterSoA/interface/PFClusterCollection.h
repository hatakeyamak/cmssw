#ifndef CUDADataFormats_PFClusterSoA_interface_PFClusterCollection_h
#define CUDADataFormats_PFClusterSoA_interface_PFClusterCollection_h

#include <vector>

#include "CUDADataFormats/PFCommon/interface/Common.h"
#include "HeterogeneousCore/CUDAUtilities/interface/HostAllocator.h"

namespace hcal {
  template <typename StoragePolicy>
  struct PFClusterCollection : public ::pf::common::AddSize<typename StoragePolicy::TagType> {
    PFClusterCollection() = default;
    PFClusterCollection(const PFClusterCollection&) = default;
    PFClusterCollection& operator=(const PFClusterCollection&) = default;

    PFClusterCollection(PFClusterCollection&&) = default;
    PFClusterCollection& operator=(PFClusterCollection&&) = default;

    typename StoragePolicy::template StorageSelector<int>::type pfc_depth;
    typename StoragePolicy::template StorageSelector<int>::type pfc_layer;
    typename StoragePolicy::template StorageSelector<int>::type pfc_detId;
    //typename StoragePolicy::template StorageSelector<int>::type pfc_neighbours;
    //typename StoragePolicy::template StorageSelector<short>::type pfc_neighbourInfos;

    typename StoragePolicy::template StorageSelector<float>::type pfc_time;
    typename StoragePolicy::template StorageSelector<float>::type pfc_energy;
    typename StoragePolicy::template StorageSelector<float>::type pfc_x;
    typename StoragePolicy::template StorageSelector<float>::type pfc_y;
    typename StoragePolicy::template StorageSelector<float>::type pfc_z;

    //m_dEta, m_dPhi, m_repCorners, backPoint

    template <typename U = typename StoragePolicy::TagType>
    typename std::enable_if<std::is_same<U, ::pf::common::tags::Vec>::value, void>::type resize(size_t size) {
      pfc_depth.resize(size);
      pfc_layer.resize(size);
      pfc_detId.resize(size);
      //pfc_neighbours.resize(8 * size);
      //pfc_neighbourInfos.resize(8 * size);

      pfc_time.resize(size);
      pfc_energy.resize(size);
      pfc_x.resize(size);
      pfc_y.resize(size);
      pfc_z.resize(size);
    }
  };  // struct PFClusterCollection

}  // namespace hcal

namespace ecal {
  template <typename StoragePolicy>
  struct PFClusterCollection : public ::pf::common::AddSize<typename StoragePolicy::TagType> {
    PFClusterCollection() = default;
    PFClusterCollection(const PFClusterCollection&) = default;
    PFClusterCollection& operator=(const PFClusterCollection&) = default;

    PFClusterCollection(PFClusterCollection&&) = default;
    PFClusterCollection& operator=(PFClusterCollection&&) = default;

    typename StoragePolicy::template StorageSelector<int>::type pfc_depth;
    typename StoragePolicy::template StorageSelector<int>::type pfc_layer;
    typename StoragePolicy::template StorageSelector<int>::type pfc_detId;
    typename StoragePolicy::template StorageSelector<int>::type pfc_neighbours;
    typename StoragePolicy::template StorageSelector<short>::type pfc_neighbourInfos;

    typename StoragePolicy::template StorageSelector<float>::type pfc_time;
    typename StoragePolicy::template StorageSelector<float>::type pfc_energy;
    typename StoragePolicy::template StorageSelector<float>::type pfc_x;
    typename StoragePolicy::template StorageSelector<float>::type pfc_y;
    typename StoragePolicy::template StorageSelector<float>::type pfc_z;

    //m_dEta, m_dPhi, m_repCorners, backPoint

    template <typename U = typename StoragePolicy::TagType>
    typename std::enable_if<std::is_same<U, ::pf::common::tags::Vec>::value, void>::type resize(size_t size) {
      pfc_depth.resize(size);
      pfc_layer.resize(size);
      pfc_detId.resize(size);
      pfc_neighbours.resize(8 * size);
      pfc_neighbourInfos.resize(8 * size);

      pfc_time.resize(size);
      pfc_energy.resize(size);
      pfc_x.resize(size);
      pfc_y.resize(size);
      pfc_z.resize(size);
    }
  };  // struct PFClusterCollection

}  // namespace ecal

#endif  // CUDADataFormats_PFClusterSoA_interface_PFClusterCollection_h
