/* Copyright 2012 Geehalel (geehalel AT gmail DOT com) */
/* This file is part of the Skywatcher Protocol INDI driver.

    The Skywatcher Protocol INDI driver is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    The Skywatcher Protocol INDI driver is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with the Skywatcher Protocol INDI driver.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <wordexp.h>
#include <math.h>

/* For IDLog/Debug */
#include <indidevapi.h>

#include "pointset.h"

double PointSet::range24(double r) {
  double res = r;
  while (res<0.0) res+=24.0;
  while (res>24.0) res-=24.0;
  return res;
}


/*
void PointSet::AltAzFromRaDec(double ra, double dec, double lst, double *alt, double *az) 
{

}
*/
 /* Using haversine: http://en.wikipedia.org/wiki/Haversine_formula */
double sphere_unit_distance(double theta1, double theta2, double phi1, double phi2) {
  double sqrt_haversin_lat = sin(((phi2 - phi1) / 2) * (M_PI / 180));
  double sqrt_haversin_long = sin(((theta2 - theta1) / 2) * (M_PI / 180));
  return (2 * asin(sqrt((sqrt_haversin_lat * sqrt_haversin_lat) +
			cos(phi1 * (M_PI / 180))*cos(phi2 * (M_PI / 180))*(sqrt_haversin_long * sqrt_haversin_long))));
}

bool compelt(PointSet::Distance d1, PointSet::Distance d2) {
  return d1.value < d2.value;
}

std::set<PointSet::Distance, bool (*)(PointSet::Distance, PointSet::Distance)> *PointSet::ComputeDistances(double alt, double az, PointFilter filter) {
  std::map<HtmID, Point>::iterator it;
  std::set<Distance, bool (*)(Distance, Distance)> *distances = new std::set<Distance, bool (*)(Distance, Distance)>(compelt);
  std::set<Distance>::iterator distit;
  std::pair<std::set<Distance>::iterator, bool> ret;
  /* IDLog("Compute distances for point alt=%f az=%f\n", alt, az);*/
  for ( it=PointSetMap->begin() ; it != PointSetMap->end(); it++ ) {
    double d = sphere_unit_distance(az, (*it).second.celestialAZ, alt, (*it).second.celestialALT);
    Distance elt;
    elt.htmID=(*it).first;
    elt.value=d;
    ret=distances->insert(elt);
    /*IDLog("  Point %lld (alt=%f az=%f): distance %f \n",  elt.htmID, (*it).second.celestialALT, (*it).second.celestialAZ, elt.value);*/
  }
  /*
  IDLog("  Ordered distances for point alt=%f az=%f\n", alt, az);
  for ( distit=distances->begin() ; distit != distances->end(); distit++ ) {
    IDLog("  Point %lld: distance %f \n",  distit->htmID, distit->value);
  }
  */
  return distances;
}

void PointSet::AddPoint(AlignData aligndata) 
{

  Point point;
  point.aligndata = aligndata;
  point.celestialAZ = (range24(point.aligndata.lst - point.aligndata.targetRA - 12.0) * 360.0) / 24.0;
  point.telescopeAZ = (range24(point.aligndata.lst - point.aligndata.telescopeRA - 12.0) * 360.0) / 24.0;
  point.celestialALT = point.aligndata.targetDEC + lat;
  point.telescopeALT = point.aligndata.telescopeDEC + lat;
  point.htmID=cc_radec2ID(point.celestialAZ, point.aligndata.targetDEC, 19);
  cc_ID2name(point.htmname,  point.htmID);
  PointSetMap->insert(std::pair<HtmID, Point>(point.htmID, point));
  /*IDLog("Adding sync point id = %lld name = %s\n ", point.htmID, point.htmname);*/
}

PointSet::Point *PointSet::getPoint(HtmID htmid) {
  return &(PointSetMap->find(htmid)->second);
}

void PointSet::Init()
{
  PointSetMap=NULL;
  PointSetXmlRoot=NULL;
}

void PointSet::Reset()
{
  if (PointSetMap) {
    PointSetMap->clear();
    delete(PointSetMap);
  }
  PointSetMap=NULL;
  if (PointSetXmlRoot)
    delXMLEle(PointSetXmlRoot);
  PointSetXmlRoot=NULL;
}

char *PointSet::LoadDataFile(const char *filename)
{
  wordexp_t wexp;
  FILE *fp;
  LilXML *lp;
  static char errmsg[512];
  AlignData aligndata;
  XMLEle *alignxml, *sitexml;
  XMLAtt *ap;
  char *sitename;
  std::map<HtmID, Point>::iterator it;
  

  if (wordexp(filename, &wexp, 0)) {
    wordfree(&wexp);
    return (char *)("Badly formed filename");
  }
  //if (filename == NULL) return;
  if (!(fp=fopen(wexp.we_wordv[0], "r"))) {
    wordfree(&wexp);
    return strerror(errno);
  }
  wordfree(&wexp);
  lp = newLilXML();
  PointSetXmlRoot = readXMLFile(fp, lp, errmsg);
  delLilXML(lp);
  if (!PointSetXmlRoot) return errmsg;
  if (!strcmp(tagXMLEle(nextXMLEle(PointSetXmlRoot, 1)), "aligndata")) return (char *)("Not an alignement data file");
  sitexml=findXMLEle(PointSetXmlRoot, "site");
  if (!sitexml) return (char *)"No site found";
  ap = findXMLAtt(sitexml, "name");
  if (ap) sitename = valuXMLAtt(ap); 
  else sitename=(char *)"No sitename";
  ap = findXMLAtt(sitexml, "lat");
  if (!ap) return (char *) "No latitude data found";
  else sscanf(valuXMLAtt(ap), "%lf", &lat);
  ap = findXMLAtt(sitexml, "lon");
  if (!ap) return (char *) "No longitude data found";
  else sscanf(valuXMLAtt(ap), "%lf", &lon);
  ap = findXMLAtt(sitexml, "alt");
  if (!ap) alt =0.0; 
  else sscanf(valuXMLAtt(ap), "%lf", &alt);
  IDLog("Align Data for site %s (lon %f lat %f alt %f)\n", sitename, lon, lat, alt);
  IDLog("  number of points: %d\n", nXMLEle(sitexml));
  PointSetMap = new std::map<HtmID, Point>();
  alignxml=nextXMLEle(sitexml, 1);
  while (alignxml) {
    //IDLog("synctime %s\n", pcdataXMLEle(findXMLEle(alignxml, "synctime")));
    sscanf(pcdataXMLEle(findXMLEle(alignxml, "synctime")), " %lf ", &aligndata.lst);
    sscanf(pcdataXMLEle(findXMLEle(alignxml, "celestialra")), "%lf", &aligndata.targetRA);
    sscanf(pcdataXMLEle(findXMLEle(alignxml, "celestialde")), "%lf", &aligndata.targetDEC);
    sscanf(pcdataXMLEle(findXMLEle(alignxml, "telescopera")), "%lf", &aligndata.telescopeRA);
    sscanf(pcdataXMLEle(findXMLEle(alignxml, "telescopede")), "%lf", &aligndata.telescopeDEC);
    //IDLog("Load alignment point: %f %f %f %f %f\n", aligndata.lst, aligndata.targetRA, aligndata.targetDEC, 
    //  aligndata.telescopeRA, aligndata.telescopeDEC);
    AddPoint(aligndata);
    alignxml=nextXMLEle(sitexml, 0);
  }
  IDLog("Resulting Alignment map;\n");
  for ( it=PointSetMap->begin() ; it != PointSetMap->end(); it++ )
    IDLog("  Point htmID= %lld: telescope alt = %f az = %f\n",  (*it).first, (*it).second.telescopeALT,  (*it).second.telescopeAZ);

  return NULL;
}
