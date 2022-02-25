#include <iostream>
using namespace std;
#include <fstream>
#include <sstream>
#include <queue>
#include <map>
#include <limits>
#include <algorithm>
#include <vector>

void convertDOT(string dfInFile, string dfOutFile) {
	// Open the file
	ifstream dfIn;
	dfIn.open(dfInFile, std::fstream::in);
	if (dfIn.fail()) throw (string("fromDOT: unable to open input file ") + dfInFile);

	ofstream dfOut;
	dfOut.open(dfOutFile, std::fstream::out);
	if (dfOut.fail()) throw (string("toDOT: unable to open output file ") + dfOutFile);

	// Print the dot to the file stream
	dfOut << "digraph " << " {" << endl;

	// Read the dot from the file stream.
	// This is a rather ad-hoc parser for a very small subset of the dot language,
	// essentially what is written by toDOT.
	// [TODO] look for a more formal parser for the dot language
	// For example: https://github.com/nkkav/gvizparse

	// A string variable used for input
	string s = "";
	string line ="";
	string op_code = "";
	string node_id;
	int nodeCount = 0;
	string edge_id;
	bool readNode = false;
	bool readEdge = false;
	map<string, int> nodeMapper;

	// Read the first line of the file; it contains the header
	std::getline(dfIn, line);
	stringstream s_stream(line);

	s_stream >> s;
	if (s != "digraph") throw (string("fromDOT: invalid graph type for this graph"));

	s_stream >> s;

	if (s != "{") { // if a name is present, discard it
		s_stream >> s;
	}

	if (s_stream.fail()) throw (string("fromDOT: syntax error, expected \"{\" after graph name"));
	if (s != "{") throw (string("fromDOT: syntax error, expected \"{\" after graph name, found " + s));

	// Now we read rest of the lines of the file
	std::getline(dfIn, line);

	while (!dfIn.eof()) {
		stringstream s_stream(line);
		
		// Check if the graph definition is done
		if (s_stream.peek() == '}') break;
		
		// Read a node id
		s_stream >> node_id;
		if (s_stream.eof()) {
			// This must be a blank line, read the next line
			getline(dfIn, line);
			// Strip all the white spaces
			//line.erase(std::remove_if(line.begin(), line.end(), ::isspace), line.end());
			continue;
		}
		if (s_stream.fail()) throw (string("fromDOT: failed to get node id"));

		// Check the next character to see if this is a node or an edge
		s_stream >> s;
		if (s == "[label") readNode = true;
		else {
			readNode = false;
			if (s == "->") readEdge = true;
			else if(s.find("[fontcolor") != string::npos) { //found fontcolor tag, skip this line
				//cout << "Skipping fontcolor attribute " << endl;
				getline(dfIn, line);
				continue;
			} else {
				readEdge = false;
				throw (string("fromDOT: syntax error, expected a node or an edge"));
			}
		}

		if (readNode) {
			nodeMapper[node_id] = nodeCount;
			nodeCount++; //for the next one
			// Read = sign
			s_stream >> s;

			//Read op
			s_stream >> s;
			
			op_code = s;
			dfOut << "   " << nodeMapper[node_id] << " [label=\"" << op_code << "\"];" << endl;
			// Add node to the graph
			//g.addNode(node_id, op_code);
		}

		if (readEdge) {
			string src_id = node_id;
			string dest_id;
			s_stream >> dest_id;
			if (s_stream.fail()) throw (string("fromDOT: failed to get destination node id"));
			dfOut << "   " << nodeMapper[src_id] << "->" << nodeMapper[dest_id] << " [label=\"1:1\"];" << endl;

		}

		// Read the next line
		getline(dfIn, line);
		// Strip all the white spaces
		//line.erase(std::remove_if(line.begin(), line.end(), ::isspace), line.end());
	}

	dfOut << "}" << endl;

	// Done, close the file
	dfOut.close();
	dfIn.close();

	// If the graph is not a dag, raise an exception
	//if ( !isDAG() ) throw (string("fromDOT: the graph is not a dag!"));
}

int main(int argc, char **argv) {
	if(argc != 3) {
		cout << "Need 2 args " << endl;
		return -1;
	}
	try {
		convertDOT(argv[1], argv[2]);
	}catch(std::string ex) {
		cout << ex << endl;
	}
	return 0;
}
