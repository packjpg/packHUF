/* -----------------------------------------------
	struct declarations
	----------------------------------------------- */

struct canonical_table {
	unsigned char* data;
	unsigned char* len_ptr;
	unsigned char* val_ptr;
	int total_size;
	int max_clen;
	int n_symbols;
	unsigned int* counts;
};

struct huffman_code {
	unsigned char len;
	unsigned int code;
};

struct huffman_conv {
	unsigned char len;
	int val;
	void* ext;
};

struct huffman_conv_set {
	huffman_conv** h;
	int hwidth;
};



/* -----------------------------------------------
	canonical table functions
	----------------------------------------------- */
	
canonical_table* build_ctable( iostream* stream );
canonical_table* build_ctable( unsigned int* counts );
canonical_table* read_ctable( iostream* stream );
void cleanup_ctable( canonical_table* ctable );


/* -----------------------------------------------
	helper functions
	----------------------------------------------- */

huffman_code** convert_ctable2hcode( canonical_table* ctable );
huffman_conv_set* convert_ctable2hconv( canonical_table* ctable );
void cleanup_hcode( huffman_code** table );
void cleanup_hconv( huffman_conv_set* conv_set );


/* -----------------------------------------------
	class for huffman decoding/data reading
	----------------------------------------------- */
	
class huffman_reader
{
	public:
	huffman_reader( unsigned char* data, int size );
	~huffman_reader( void );
	// public functions
	int decode_symbol( huffman_conv_set* table );
	unsigned int read_bits( int n );
	// not the fastest method to read single bits/bytes
	// luckily, we won't need it that often
	unsigned char read_bit( void );
	int getpos();
	
	private:
	// utility functions
	inline void advance_bitstream( int n );
	inline void advance_bitstream_1( void );
	
	// storage
	abitreader* bit_reader;
	unsigned int bit_buffer;
};


/* -----------------------------------------------
	class for huffman encoding/data writer
	----------------------------------------------- */
	
class huffman_writer
{
	public:
	huffman_writer( int adds );
	~huffman_writer( void );
	// public functions
	void encode_symbol( huffman_code** hcodes, int symbol );
	void write_bits( unsigned int bits, int n );
	void write_bit( unsigned char bit );
	unsigned char* getptr( void );
	int getpos();
	
	private:
	// storage
	abitwriter* bit_writer;
};
