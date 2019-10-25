/*******************************************************************************
    PROFASI: Protein Folding and Aggregation Simulator, Version 1.5
    Copyright (C) (2012)  Anders Irback and Sandipan Mohanty
    Email: profasi@thep.lu.se
    Home Page: http://cbbp.thep.lu.se/activities/profasi/
    Version control (git) : https://trac.version.fz-juelich.de/PROFASI

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License
    (see PROFASI/gpl.txt).

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.
********************************************************************************/

#ifndef PRF_XML_HH
#define PRF_XML_HH
#include <string>
#include <map>
#include <deque>
#include <list>
#include "profasi_io.hh"

/**
* \defgroup profasi_XML ProFASi XML module
* \ingroup utilities
* @brief Basic XML parsing for ProFASi
*
* The namespace prf_xml is a collection of classes and functions meant to
* provide basic support for parsing XML in ProFASi. It is not meant to be
* a comprehensive XML parser module. There are more established tools for
* doing such things. The purpose of this module is to provide the minimal
* functionality that is required in ProFASi while still keeping the package
* self contained.
*
* One non-trivial feature of this module is its ability to autogenerate
* XML nodes using a given template. See the documentation of
* \ref prf_xml::XML_Node::interpret_formatted_data() for details.
*
* In this preliminary version, special XML tags such as \<?xml ...?\> and
* DOCTYPE declarations are collected, but then ignored. In the future the
* module may be extended, upon necessity, to process such tags.
*
* An XML tree is made of inter-related XML nodes. Each node has a name, some
* "attributes", and possibly, one or more "child nodes". A tree has to contain
* one and only one "root" node. Every other node in the tree can be reached
* by starting from the root node and recursively browsing the child nodes.
*/
namespace prf_xml
{
    //forward declarations
    class Chunk;
    class XML_Node;
    typedef enum {ERROR, OK} state_type;
    typedef enum {bad, text, begin_tag, end_tag, special_tag} chunk_type;
    //! A function returning an XML tree taking an XML file name
    /**
    * This is probably how this module will most often be used. You have an
    * XML file. Pass it to this function and get a pointer to the root
    * node of the XML tree in that file. A lot of things happen in between,
    * such as syntax checking of the XML tags, but such things are normally
    * not of interest to the outside world. Either the file parses nicely
    * into an XML tree or it does not. In case of errors, the returned pointer
    * will be NULL.
    *
    * A little note on memory management: this function allocates memory as
    * needed for the root and the child XML nodes. It is the responsibility
    * of the outside program to destroy the newly allocated node by calling
    * \c delete on it.
    * \sa XML_Node
    * \ingroup profasi_XML
    */
    XML_Node *get_xml_tree(std::string xml_file_name);

    //! Convert a string representation of an XML tree into a tree
    XML_Node *make_xml_tree(std::string xml_node_text);

    //! A class to represent an XML node
    /**
    * An XML node represents the information enclosed within the scope of
    * an XML tag. It has a name, some data, and a deque of pointers to other
    * nodes regarded as its children. These child nodes are the XML tags
    * contained within the scope of the original tag. For example, \n \n
    *
    \verbatim
    <protein>
    <sequence>
    ACE * KLVFFAE * NH2
    </sequence>
    <group index="0" type="ACE">
    <dof>
    4.0092316407542254
    </dof>
    </group>
    <group index="1" type="LYS">
    <dof>
    -0.3131207777043232   -1.3363358975000719   -3.1415926535897931
    4.5758682133671185   2.8279753720956338   1.1606113829381837
    3.1121544489086883   1.9936435953881690
    </dof>
    </group>
    </protein>
    \endverbatim
    *
    * represents a short tree with information about the beginning part of a
    * peptide. The root node in this case has a name "protein". It has a few
    * child nodes (i) A child node called "sequence", which contains the
    * sequence of the peptide as its data and no children of its own. (ii) 2
    * child nodes with name "group". In this instance, these two are
    * distinguishable by their attributes "index" and "type". Each group
    * contains a child node of its own, called "dof", with the double precision
    * values of, presumably, the degrees of freedom in that residue. Notice
    * that the name of a node appears in the opening and closing tags that
    * define the scope of the node. Tags for the child nodes are enclosed
    * within the scope of the parent node. The opening tag can optionally
    * contain one or more attributes of arbitray names. Attribute values
    * must be enclosed in quotes. There are no commas between different
    * attributes, so that they are only distinguished by identifying
    * name="value" patterns.
    *
    * This class is an interface to the attributes and information contained
    * in an XML node and its children.
    * \ingroup profasi_XML
    */
    class XML_Node
    {
    public:
        //! Default constructor
        XML_Node();
        ~XML_Node();
        //! Copy constructor
        XML_Node(const XML_Node&);
        //! Create a real detached copy with the same information
        /**
          The copy constructor copies all elements, including the pointers.
          The child nodes are pointers to other XML_Nodes. So, if we construct
          a copy using the copy constructor, we end up with two XML_Nodes which
          have elements pointing to the same memory address. This could lead
          to problems when the nodes are deleted. When we want a duplicate,
          one should use this function instead.
          */
        void make_clone_of(const XML_Node *);
        //! Construct from a range of "Chunks"
        /**
        * Normally this should not concern you unless you intend to want to
        * improve the construction of the XML tree in some way. It takes a
        * range of Chunk objects and tries to generate a node out of it. The
        * range of chunks is assumed to start with a begin tag, end with an
        * end tag and contain only data within these tags, i.e., no enclosed
        * child nodes. Read the documentation of the Chunk and XML_Mini
        * classes for further details of how and why it is used. If you only
        * want to use the XML module, you can safely ignore this function.
        *
        * \sa Chunk
        * \sa XML_Mini
        */
        XML_Node(std::list<prf_xml::Chunk>::iterator bg,
                 std::list<prf_xml::Chunk>::iterator nd);
        //! Assignment operator
        XML_Node &operator=(const XML_Node &);
        //! Create a node with a tag name and some data
        /**
        * @param nm Name or tag name for the XML node
        * @param vl Text data for the node
        */
        XML_Node(std::string nm, std::string vl);
        //! Set the tag name (no syntax check!!)
        inline void set_name(std::string nm) {the_name=nm; }
        //! Set the text data (no syntax check!!)
        inline void set_value(std::string vl) { the_value=vl; }
        inline void set_attribute(std::string atrnm, std::string atrvl) {
            attr[atrnm]=atrvl;
        }
        //! Retrieve the name of a node
        inline std::string name() const { return the_name; }
        //! Retrieve the text data stored in the node
        /**
        * If data is stored in the file non-contiguously with child XML tags
        * interspersed between text data blocks, the data belonging to the
        * node (and not to the child nodes) will be collected in one place.
        * This function then returns all the data belonging to the node.
        */
        inline std::string value() const { return the_value; }
        //! Whether the node is a leaf node
        /**
          Returns true, if the node has no children, no attributes, and the
          value text does not contain more than 80 characters.
          */
        inline bool is_leaf_node() const {
            return chldrn.empty() && childptr.empty() && the_value.size()<80;
        }
        //! Retrieve a particular attribute
        /**
        * Empty string if the desired attribute does not exist
        */
        inline std::string attribute(std::string a) { return attr[a]; }
        //! Number of attributes
        inline size_t n_attributes() const { return attr.size(); }
        //! Reference to the attribute map
        inline std::map<std::string,std::string> &attributes() { return attr; }
        //! Number of children
        inline size_t n_children() const { return chldrn.size(); }
        //! i'th child node
        inline XML_Node * child(size_t i) { return chldrn[i]; }
        //! Pointer to the first child node with a given name
        XML_Node * child(std::string nm);
        //! Reference to the deque of children
        inline std::deque<XML_Node *> &children() {return chldrn;}
        //! Add a child node
        int add_child_node(XML_Node * ch);
        //! Remove a node from the deque of child nodes
        int remove_child_node(std::string nm);
        //! Remove a child node
        int remove_child_node(XML_Node *nd);
        //! Disown a child node
        /**
          * If nd is a child node, it will be removed from the child list but
          * its contents will not be touched. So, the responsibility of deleting
          * that node will no longer be with the original parent node.
          */
        int disown_child_node(XML_Node *nd);
        //! Add a new child node with given name and value
        int add_child_node(std::string nm, std::string vl);
        //! Clear all child nodes
        void clear_child_nodes();
        //! Whether or not the XML node is in an error free state
        inline state_type state() const { return mystate; }
        //! Whether it encloses another XML node
        inline bool encloses(XML_Node &gn) {
            return cmin<=gn.cmin and cmax>=gn.cmax;
        }
        //! Make a string out of the information in the node and its children
        std::string make_string() const ;
        //! Write the XML node and its children for visualisation
        /**
        * This function is only to see that the XML parsing has worked. It
        * writes the node tag, defines its scope like the scope of a C function,
        * shows the attributes like arguments to the function, uses indentation
        * to represent the hierarchy.
        *
        */
        void write(prf::Output &op, int indent_level=0);
        //! Interpret formatted data as a series of child nodes
        /**
          This function interprets formatted data, and hence implements the
          XML template handling for ProFASi.
          \sa xml_templates
          */
        int interpret_formatted_data();
    private:
        std::string the_name,the_value;
        std::map<std::string, std::string> attr;
        std::multimap<std::string, XML_Node *> childptr;
        std::deque<XML_Node *> chldrn;
        state_type mystate;
        size_t cmin, cmax; //For internal use only
    };

    //! A "chunk" of data in an XML file
    /**
    * A chunk is a string of characters with the following properties:\n \n
    * \li Either it starts with a '&lt;' and ends with a '&gt;' character
    * \li or it does not contain '&lt;' or '&gt;' anywhere.
    *
    * An XML file is a sequence of chunks. Depending on its content a chunk
    * of data can be classified as an XML begin tag, an end tag, text data,
    * a special tag, or garbage that can not be parsed as proper XML. This
    * class provides ways to handle chunks. It is a helper class for
    * constructing XML node trees out of an XML file.
    *
    * A Chunk object not only stores the text it represents, it can keep track
    * of its location in the file. For instance, if a chunk &lt;group&gt;
    * appears as characters 123 through 129 in a file, it should be told that
    * its limits are 123 to 130. When XML nodes are constructed from a range
    * of chunks, these limits come in handy. They define the scope of one
    * XML node. Using that scope, it is possible to tell if an XML node is a
    * child of another.
    * \sa XML_Node, XML_Mini
    * \ingroup profasi_XML
    */
    class Chunk
    {
    public:
        //! Default constructor
        Chunk();
        ~Chunk();
        //! Copy constructor
        Chunk(const Chunk &);
        //! Construct using text and limits
        Chunk(std::string txt, size_t bg, size_t nd);
        //! Assignment operator
        Chunk &operator=(const Chunk &);
        //! Where it starts
        inline size_t begin() const { return i1; }
        //! Where it ends
        inline size_t end() const { return i2; }
        //! A line number, if it was derived from a file
        inline size_t line_number() const { return line_no; }
        //! Assign a line number
        inline void line_number(size_t i) {line_no=i;}
        //! What kind of chunk it is
        /**
        * All chunks are created as type text. This function does not
        * analyze the chunk data to see what type it is. It only reports
        * the stored value for its type. The function determine_type() must
        * be called before using the function type().
        */
        inline chunk_type type() const { return mytype; }
        //! Whether it is within a given scope
        inline bool is_in(size_t r1, size_t r2) const {
            return (i1>=r1 && i1<r2 && i2<=r1 &&i2<r2);
        }
        //! Set the chunk text
        inline void set_text(std::string s) { the_text=s; }
        //! Set the limits
        inline void set_limits(size_t i, size_t j) { i1=i; i2=j; }
        //! The text data of the chunk
        inline std::string show() { return the_text; }
        //! Info about type, limits and text data
        std::string info();
        //! Determine what kind of a chunk it is
        /**
        * Determines whether the chunk is an XML begin tag, an end tag, a
        * special tag, data text or garbage. It uses helper functions
        * syntax_check_text() to validate data text, valid_name() to validate
        * the tag name of the chunk or the attributes and
        * syntax_check_attributes() to determine if the attributes have been
        * specified with proper syntax. A valid end tag is of the form
        * &lt;/a_valid_name&gt;, where the validity of a_valid_name is
        * determined using the function valid_name(). Spaces between the end
        * of a_valid_name and the closing '&gt;' character are tolerated. A
        * valid begin tag is of the form
        * &lt;a_valid_name valid_attribute_list&gt;.
        * \sa valid_name(), syntax_check_text(), syntax_check_attributes()
        */
        chunk_type determine_type();
        //! Return the tag name (name sans the &lt; &gt; signs)
        std::string tag_name();
        //! If it is a begin tag, return the specified attributes
        std::map<std::string,std::string> attribute_list();
    private:
        std::string the_text, tagname;
        size_t i1,i2,line_no;
        chunk_type mytype;
        std::map<std::string,std::string> atl;
    };

    //! A small class to do the book keeping for XML parsing
    /**
    * Divides the task into simple logical blocks like reading the file into
    * a list of chunks and then building the XML tree from the chunks. Not
    * really intended for extensive use. Use the interface function
    * get_xml_tree() instead. This is a helper class for the implementation
    * of the functionality in get_xml_tree()
    * \sa get_xml_tree()
    * \ingroup profasi_XML
    */
    class XML_Mini
    {
    public:
        //! Default constructor
        XML_Mini();
        ~XML_Mini();
        //! Read XML file as a list of Chunk s
        /**
        * Reads an XML file specified by its name and makes a list of Chunk
        * objects out of it. Each of these chunks, by construction, either does
        * not contain any &lt; or &gt; tags, or they appear simultaneously as
        * the first and last characters respectively.
        */
        size_t read_file(std::string fl);
        //! Build XML tree and return pointer to the root element
        /**
        * Constructs an XML tree out of a list of chunks perviously collected
        * using read_file(). Here is the algorithm:
        *
        * <ol>
        * <li> Check that no chunk contains errors
        * <li> Remove special tags to a separate list
        * <li> Initialize empty node list "node_buffer"
        * <li> Search for sequences of tags that read begin, text ... end. In
        * other words, a begin followed by an end with only text in between.
        * <li> Construct a node using this range of chunks, and remove them from
        * the list.
        * <li> If the constructed node contains errors (begin and end
        * tag names don't match, name or attribute syntax invalid or the
        * enclosed text has errors) delete the newly created node and put
        * parser state to ERROR.
        * <li> If an error free node was successfully created, browse the
        * node_buffer for previously created nodes to find nodes enclosed by
        * the newly created node. If any are found, add them as children and
        * remove them from the node_buffer. Finally add the new node to the
        * node_buffer
        * <li> Repeat everything from step 4 onwards. Go back to the start
        * of the survibing chunk list if the end of the list is encountered and
        * the parser is not in an ERROR state. Continue until the list of
        * Chunks is empty or the process stagnates.
        * <li> If there is exactly one element in the node_buffer, return it as
        * the root element. Else return NULL
        * </ol>
        */
        XML_Node *build_tree();
        //! Takes a string, and parses it into a list of Chunk objects
        size_t parse_data(std::string &dat);
        //! Clear all data
        void clear();
    public:
        std::list<Chunk> data;
    };
    //! Check if a given string is a valid name for a tag
    /**
    * A valid name contains only alpha-numeric characters or the underscore '_'.
    * It can not start with a digit or the underscore.
    */
    bool valid_name(std::string nm);
    //! Check if a given string is valid text for XML
    /**
    * Check that the text does not contain any &lt; or &gt; characters. Although
    * strictly speaking &gt; is not forbidden in the XML specifications, we
    * do not allow it for simplicity. Both these characters must be written
    * using entity references described below. The character &amp; should only
    * appear in the following combinations: &amp;lt;, &amp;gt;, &amp;amp;,
    * &amp;qot; and &amp;apos;. These expressions are called entity references
    * representing the characters &lt; &gt; &amp; &qot; and &apos; respectively.
    *
    */
    bool syntax_check_text(std::string txt);
    //! Check if a given string is a valid list of attributes
    /**
    * To be a good attribute list, a string must be a white space separated
    * list of  <i>attribute_name = "attribute_value"</i> blocks. The
    * attribute names must be valid names as defined under valid_name().
    * Attribute values must be quoted, and must obey the syntax rules for
    * good XML text.
    */
    bool syntax_check_attributes(std::string att,
                                 std::map<std::string, std::string> &atl);
    bool get_atr_name(std::string att,size_t &ipos, std::string &nm);
    bool get_atr_text(std::string att, size_t &ipos, std::string &tx);
    /**
  \page xml_templates XML templates: Auto-generate complex XML nodes from tabular data

  Often XML data consists of a large sequence of nodes of identical
  format with different entries. Although it is more convenient to
  analyze the data as a sequence of XML nodes, writing out each
  node in the standard XML format is unnecessarily verbose. Storing
  such data in the table form in a file with comma or white space
  separated columns is space efficient and convenient from the
  point of view of scripting. But this approach makes us commit to
  a certain data layout and not very accommodative of unexpected
  new directions. ProFASi's XML module provides a template handling
  mechanism to auto-generate XML nodes from tabular data.

  An XML_Node with a child node of name "formatted_data" can
  substitute a series of child nodes based on the contents of a
  \e formatted_data child. Each \e formatted_data node must have
  children of names "format" and "data" (See example below). Each
  line of input inside the \e data field will be interpreted as
  a record, represented by XML nodes described by the \e format.
  The format field can address different fields in the data line
  in \e awk style, \e $1 \e $2 etc. Example:

  \verbatim
  <some_node>
  <arbitrary_property_1>
    This is some generic property.
  </arbitrary_property_1>
  <formatted_data>
    <format name="snapshot">
      <MC_time>$1</MC_time>
      <temperature>$2</temperature>
      <energy>$3</energy>
      <helix>$4</helix>
      <strand>$5</strand>
    </format>
    <data>
      999  0  33.391879  0.071429  0.357143
      1999  3  67.611951  0.000000  0.071429
      2999  2  53.211118  0.000000  0.000000
    </data>
  </formatted_data>
  </some_node>
  \endverbatim

  If function prf_xml::XML_Node::interpret_formatted_data() is called on XML_Node
  \tt some_node, the content of the node \tt some_node will change
  to the following:
  \verbatim
  <some_node>
  <arbitrary_property_1>
    This is some generic property.
  </arbitrary_property_1>
  <snapshot>
      <MC_time>999</MC_time>
      <temperature>0</temperature>
      <energy>33.391879</energy>
      <helix>0.071429</helix>
      <strand>0.357143</strand>
  </snapshot>
  <snapshot>
      <MC_time>1999</MC_time>
      <temperature>3</temperature>
      <energy>67.611951</energy>
      <helix>0.000000</helix>
      <strand>0.071429</strand>
  </snapshot>
  <snapshot>
      <MC_time>2999</MC_time>
      <temperature>2</temperature>
      <energy>53.211118</energy>
      <helix>0.000000</helix>
      <strand>0.000000</strand>
  </snapshot>
  </some_node>
  \endverbatim
  The advantage of being able to do this is that an application can
  process the data by accessing XML nodes with given names. The
  correspondence between those nodes and columns of data in a
  tabular file is left completely open. So, an external application
  can generate data in arbitrary format. So long as there is a
  format specifier block added and the data copied into a data
  field, it can be interpreted by the processing program without
  any change.

  Tip: In the formatted data node, there is an alternative
  to the data node, called \tt import_data. The content of
  \tt import_data node should be a file name. During interpretation,
  the block
  \verbatim
  <import_data>tabular_file.dat</import_data>
  \endverbatim
  is interpreted as if it was
  \verbatim
  <data>
  contents of the file tabular_file.dat
  </data>
  \endverbatim
  This way, you can construct many small XML files of arbitrarily
  different node structure using the information in the same
  tabular data file.

  The interpretation of formatted data happens only if one calls the function
  interpret_formatted_data() for the specific XML node. This happens by default
  for all ProFASi simulation programs which read XML configurations. But this
  happens because they explicitly call the interpret_... function after
  retrieving an XML node from the XML file. If you want to use this feature
  in a new program you write, make sure you do something like this:
  \verbatim
  prf_xml::XML_Node * root=prf_xml::get_xml_tree("some_file.xml");
  root->interpret_formatted_data();
  \endverbatim

  */
}

#endif
