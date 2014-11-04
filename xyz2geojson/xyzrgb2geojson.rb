#!/usr/bin/env ruby

# Copyright (c) 2011-2013, James Klassen
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.

def header()
	print <<EOF
{ "type": "FeatureCollection", 
        "features": [
EOF
end
	
def feature( x,y,z,r,g,b,i )
	print <<EOF
	{ 
                "type":"Feature",
                "geometry": {   
                        "type": "Point", 
                        "coordinates": [#{x},#{y},#{z}]
                },
                "properties": {
                        "X":"#{x}",
                        "Y":"#{y}",
                        "Z":"#{z}",
                        "Intensity":"#{i.to_i.to_s}",
                        "Classification":"2",
                        "PointSourceId":"20",
                        "Red":"#{r.to_i.to_s}",
                        "Green":"#{g.to_i.to_s}",
                        "Blue":"#{b.to_i.to_s}"
                }
        }
EOF
end

def footer()
	print <<EOF
	]
}
EOF
end


def read_records() 
	fields = nil
	STDIN.each_line do |line|
		line.chomp!

		if(fields)
			yield(fields)
			print ",\n"
		end
	
		fields = line.split(" ")
	end
	yield(fields) if fields
end


begin
	header
	read_records do |fields|
		feature(*fields)
	end	
	footer
end
