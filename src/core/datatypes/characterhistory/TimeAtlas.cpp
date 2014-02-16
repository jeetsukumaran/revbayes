//
//  TimeAtlas.cpp
//  rb_mlandis
//
//  Created by Michael Landis on 12/3/13.
//  Copyright (c) 2013 Michael Landis. All rights reserved.
//

#include "TimeAtlas.h"
#include "TimeAtlasDataReader.h"

using namespace RevBayesCore;

TimeAtlas::TimeAtlas(TimeAtlasDataReader* tadr) : areas(tadr->getAreas()), times(tadr->getTimes())
{

}

TimeAtlas::TimeAtlas(const TimeAtlas& a)
{
    *this = a;
}

TimeAtlas& TimeAtlas::operator=(const TimeAtlas& a)
{
    if (this != &a)
    {
        times = a.times;
        for (size_t i = 0; i < a.areas.size(); i++)
        {
            std::vector<GeographicArea*> tmpAreas;
            for (size_t j = 0; j < a.areas[i].size(); j++)
            {
                tmpAreas.push_back(new GeographicArea(*a.areas[i][j]));
            }
            areas.push_back(tmpAreas);
        }
    }
    
    return *this;
}

TimeAtlas* TimeAtlas::clone(void) const
{
    return new TimeAtlas(*this);
}

std::vector<double> TimeAtlas::getTimes(void)
{
    return times;
}

std::vector<std::vector<GeographicArea*> > TimeAtlas::getAreas(void)
{
    return areas;
}