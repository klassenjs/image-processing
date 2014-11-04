/*
 * Copyright (c) 2011-2013, James Klassen
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 *
 * Change History:
 *   10/25/2012 - modify xyz2pcd to grab color information from image3 output
 *   10/26/2012 - reorg for one pass through input file
 */

#include <iostream>
#include <pcl/io/pcd_io.h>
#include <pcl/point_types.h>
#include <stdio.h>

int
  main (int argc, char** argv)
{
  pcl::PointCloud<pcl::PointXYZRGB> cloud;

  // Fill in the cloud data
  cloud.width    = 1000;
  cloud.height   = 1;
  cloud.is_dense = false;
  cloud.points.resize (cloud.width * cloud.height);

  FILE* xyz = fopen(argv[1], "r");

  size_t i = 0;
  float x,y,z,r,g,b,ir;
  while( !feof(xyz) )
  {
    fscanf( xyz, "%f %f %f %f %f %f %f\n", &x, &y, &z, &r, &g, &b, &ir );
    if( i >= cloud.width ) {
	cloud.width += 1000;
	cloud.points.resize (cloud.width * cloud.height);
    }
    cloud.points[i].x = x;
    cloud.points[i].y = y;
    cloud.points[i].z = z;

    cloud.points[i].r = r;
    cloud.points[i].g = g;
    cloud.points[i].b = b;
    i++;
  }
  fclose(xyz);
  cloud.width = i;
  cloud.points.resize (cloud.width * cloud.height );

  pcl::io::savePCDFileBinary ("test_pcd.pcd", cloud);
  std::cerr << "Saved " << cloud.points.size () << " data points to test_pcd.pcd." << std::endl;

  return (0);
}
