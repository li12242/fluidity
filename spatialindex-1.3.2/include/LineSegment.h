// Spatial Index Library
//
// Copyright (C) 2004  Navel Ltd.
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
//
//  Email:
//    mhadji@gmail.com

#ifndef __spatialindex_linesegment_h
#define __spatialindex_linesegment_h

namespace SpatialIndex
{
	class LineSegment : public IObject, public virtual IShape
	{
	public:
		LineSegment();
		LineSegment(const double* startPoint, const double* endPoint, size_t dimension);
		LineSegment(const Point& startPoint, const Point& endPoint);
		LineSegment(const LineSegment& l);
		virtual ~LineSegment();

		virtual LineSegment& operator=(const LineSegment& p);
		virtual bool operator==(const LineSegment& p) const;

		//
		// IObject interface
		//
		virtual LineSegment* clone();

		//
		// ISerializable interface
		//
		virtual size_t getByteArraySize();
		virtual void loadFromByteArray(const byte* data);
		virtual void storeToByteArray(byte** data, size_t& length);

		//
		// IShape interface
		//
		virtual bool intersectsShape(const IShape& in) const;
		virtual bool containsShape(const IShape& in) const;
		virtual bool touchesShape(const IShape& in) const;
		virtual void getCenter(Point& out) const;
		virtual size_t getDimension() const;
		virtual void getMBR(Region& out) const;
		virtual double getArea() const;
		virtual double getMinimumDistance(const IShape& in) const;

		virtual double getMinimumDistance(const Point& p) const;
		//virtual double getMinimumDistance(const Region& r) const;
		virtual double getRelativeMinimumDistance(const Point& p) const;
		virtual double getRelativeMaximumDistance(const Region& r) const;
		virtual double getAngleOfPerpendicularRay();

		virtual void makeInfinite(size_t dimension);
		virtual void makeDimension(size_t dimension);

	public:
		size_t m_dimension;
		double* m_pStartPoint;
		double* m_pEndPoint;

		friend class Region;
		friend class Point;
		friend std::ostream& operator<<(std::ostream& os, const LineSegment& pt);
	}; // Point

	std::ostream& operator<<(std::ostream& os, const LineSegment& pt);
}

#endif /*__spatialindex_linesegment_h*/
