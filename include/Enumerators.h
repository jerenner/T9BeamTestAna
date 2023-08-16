#include <iostream>
#include <vector>
#include <unordered_map>

//Author: Junjie Xia
//First version: 19 Jul 2023 

//from: https://stackoverflow.com/questions/18837857/cant-use-enum-class-as-unordered-map-key
// pass as third argument to unordered_map when using an enum class as key (needed e.g. for gcc 4.8.5)
struct EnumClassHash
{
    template <typename T>
    std::size_t operator()(T t) const
    {
        return static_cast<std::size_t>(t);
    }
};

// define for which kind of event the target is, e.g. mupi 1pulse, etc.
enum class EventRegion
{
    not_defined = -1,
    mupi_1pulse = 0,
    e_npulses = 1,
    mupi_npulses = 2,
    zeroPeak = 3,
    onePeak = 4,
    multiPeaks = 5,
    EventRegionCount
};

std::unordered_map<EventRegion, std::string, EnumClassHash> EventRegionNames =
  {
   {EventRegion::mupi_1pulse, "mupi_1plike"},
   {EventRegion::e_npulses, "e_nplike"},
   {EventRegion::mupi_npulses, "mupi_nplike"},   
   {EventRegion::zeroPeak, "zeroPeak"},   
   {EventRegion::onePeak, "onePeak"},   
   {EventRegion::multiPeaks, "multiPeaks"},   
  };

enum class DetChannels
{
    undefined = -1,
    act_e_veto0 = 0,
    act_e_veto1,
    act_1_R = 2,
    act_1_L,
    act_2_R = 4,
    act_2_L,
    act_3_R = 6,
    act_3_L,
    tof_0_0 = 8,
    tof_0_1,
    tof_0_2,
    tof_0_3,
    tof_1_0,
    tof_1_1,
    tof_1_2,
    tof_1_3,
    hc0,
    hc1,
    pbg = 18,
    nDetChannels,
};

std::unordered_map<DetChannels, std::string, EnumClassHash> DetChNames =
{
   {DetChannels::act_e_veto0, "act0_r"},
   {DetChannels::act_e_veto1, "act0_l"},
   {DetChannels::act_1_R,     "act1_r"},
   {DetChannels::act_1_L,     "act1_l"},
   {DetChannels::act_2_R,     "act2_r"},
   {DetChannels::act_2_L,     "act2_l"},
   {DetChannels::act_3_R,     "act3_r"},
   {DetChannels::act_3_L,     "act3_l"},
   {DetChannels::tof_0_0,     "trg00"},
   {DetChannels::tof_0_1,     "trg01"},
   {DetChannels::tof_0_2,     "trg02"},
   {DetChannels::tof_0_3,     "trg03"},
   {DetChannels::tof_1_0,     "trg10"},
   {DetChannels::tof_1_1,     "trg11"},
   {DetChannels::tof_1_2,     "trg12"},
   {DetChannels::tof_1_3,     "trg13"},
   {DetChannels::hc0,         "hole0"},
   {DetChannels::hc1,         "hole1"},
   {DetChannels::pbg,         "Pbg"},
};

std::vector<std::pair<int, int>> vec_Digitizers;

std::unordered_map<DetChannels, std::pair<int, int>, EnumClassHash> DCtoDGmap;
