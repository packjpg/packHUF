#include <stdlib.h>
#include <string.h>
#include "bitops.h"
#include "huffcoder.h"

// maximum code length, must be 1 <= X < sizeof( int )
#define MAX_CLEN ( (int) ( (sizeof(int)*8) - 1 ) )
// maximum hwidth - must be 1 <= X < 31
// higher mean more memory consumption and speed
#define MAX_HWIDTH 8
// bit buffer size - careful changing this!
// should not be smaller than 8 or MAX_HWIDTH
#define BIT_BUFFER_SIZE	16


/* -----------------------------------------------
	additional declarations for ctable functions
	----------------------------------------------- */

struct huffman_node {
	unsigned int count;
	// int clen;
	int val;
	huffman_node* next[2];
};

static void quicksort_nodes( huffman_node** lbound, huffman_node** rbound );
static void count_lengths( huffman_node* node, unsigned char* clen );
static void cleanup_tree( huffman_node* node );




/* -----------------------------------------------
	build canonical huffman table
	----------------------------------------------- */	
canonical_table* build_ctable( iostream* stream ) {
	unsigned int counts[ 256 ] = { 0 };
	unsigned char v;
	
	
	// count symbol occurences
	stream->rewind();
	while ( stream->read( &v, 1, 1 ) == 1 )
		counts[ (int) v ]++;
	
	
	return build_ctable( counts );
}


/* -----------------------------------------------
	build canonical huffman table
	----------------------------------------------- */	
canonical_table* build_ctable( unsigned int* counts ) {
	unsigned char clen[ 256 + 1 ] = { 0 };
	canonical_table* ctable;
	huffman_node** nodes;
	huffman_node** tree;
	huffman_node* n;
	int nn, ml;
	int i, j;
	
	
	// build nodes
	nodes = ( huffman_node** ) calloc( 256+1, sizeof( huffman_node* ) );
	if ( nodes == NULL ) return NULL;
	nodes[ 256 ] = NULL;
	for ( i = 0; i < 256; i++ ) {
		n = ( huffman_node* ) calloc( 1, sizeof( huffman_node ) );
		if ( n == NULL ) return NULL;
		n->count = 0;
		n->val = i;
		n->next[0] = NULL;
		n->next[1] = NULL;
		n->count = counts[ i ];
		nodes[i] = n;
	}
	
	// quicksort by counts
	quicksort_nodes( nodes, nodes + 255 );
	
	// remove zero count nodes (at least one node has to be left over)
	for ( i = 255; i > 0; i-- ) {
		if ( nodes[i]->count > 0 ) break;
		free( nodes[ i ] );
		nodes[i] = NULL;
	}	
	nn = i + 1; // number of remaining nodes (must be >= 1)
	
	// make a copy of the current nodes array
	tree = ( huffman_node** ) calloc( nn + 1, sizeof( huffman_node* ) );
	if ( tree == NULL ) return NULL;
	tree[ nn ] = NULL;
	memcpy( tree, nodes, nn * sizeof( huffman_node* ) );
	
	// build huffman tree
	for ( i = nn - 1; i > 0; i-- ) {
		// build new node combining two nodes with smallest count
		n = ( huffman_node* ) calloc( 1, sizeof( huffman_node ) );
		if ( n == NULL ) return NULL;
		n->count = tree[i]->count + tree[i-1]->count;
		n->val = -1;
		n->next[0] = tree[i-1];
		n->next[1] = tree[i];
		// replace old nodes, sort new node to its place
		for ( j = i - 1; ( j > 0 ) && ( n->count > tree[j-1]->count ) ; j-- )
			tree[j] = tree[j-1];
		tree[i] = NULL;
		tree[j] = n;
	}
	
	// count code lengths using tree
	count_lengths( (*tree), clen );
	if ( clen[0] == 1 ) {
		clen[0] = 0;
		clen[1] = 1;
	}
	
	// find max code length & limit code length if needed
	for ( ml = 256; ( clen[ml] == 0 ) && ( ml > 0 ); ml-- );
	while ( ml > MAX_CLEN ) {
		while ( clen[ml] > 0 ) {
			// find end node at least 2 below
			for ( i = ml - 2; clen[i] == 0; i-- );
			// replace a tree node with a value
			clen[ml]   -= 2;
			clen[ml-1] += 1;
			// replace a value with a tree node
			clen[i+1]  += 2;
			clen[i]    -= 1;
		}
		while( clen[--ml] == 0 );
	}
	if ( ml == 0 ) ml = ( nn == 256 ) ? 8 : 1;
	
	// build canonical table
	ctable = ( canonical_table* ) calloc ( 1, sizeof( canonical_table ) );
	if ( ctable == NULL ) return NULL;
	if ( ( nn == 256 ) && ( ml == 8 ) ) { // incompressible case
		// also handle this case in huffcoder routines!
		ctable->max_clen = 8;
		ctable->n_symbols = 256;
		ctable->total_size = 2;
		ctable->len_ptr = NULL;
		ctable->val_ptr = NULL;
		ctable->data = ( unsigned char* ) calloc( 2, 1 );
		ctable->counts = ( unsigned int* ) calloc( nn, sizeof( int ) );
		if ( ( ctable->data == NULL ) || ( ctable->counts == NULL ) ) return NULL;
		for ( i = 0; i < 256; i++ ) ctable->counts[ i ] = nodes[i]->count;
		ctable->data[0] = 0xFF;
		ctable->data[1] = 0xFF;
	} else { // generic routine
		ctable->max_clen = ml;
		ctable->n_symbols = nn;
		ctable->total_size = 1 + ml + nn;
		ctable->data = ( unsigned char* ) calloc( ctable->total_size, 1 );
		ctable->counts = ( unsigned int* ) calloc( nn, sizeof( int ) );
		if ( ( ctable->data == NULL ) || ( ctable->counts == NULL ) ) return NULL;
		ctable->len_ptr = ctable->data + 1;
		ctable->val_ptr = ctable->data + 1 + ml;
		ctable->data[0] = ml;
		memcpy( ctable->len_ptr, clen + 1, ml );
		for ( i = 0; i < nn; i++ ) {
			ctable->val_ptr[ i ] = nodes[i]->val;
			ctable->counts[ i ] = nodes[i]->count;
		}
	}
	
	
	// free memory
	cleanup_tree( (*tree) );
	free( tree );
	free( nodes );
	
	
	return ctable;
}


/* -----------------------------------------------
	read canonical huffman table from file
	----------------------------------------------- */
canonical_table* read_ctable( iostream* stream ) {
	canonical_table* ctable;
	unsigned char ml;
	int nn;
	int t;
	
	
	// allocate memory
	ctable = (canonical_table*) calloc( 1, sizeof( canonical_table ) );
	if ( ctable == NULL ) return NULL;
	ctable->counts = NULL;
	
	// read first byte, allocate memory
	stream->read( &ml, 1, 1 );
	ctable->data = (unsigned char*) calloc( 1 + ml, 1 );
	if ( ctable->data == NULL ) {
		free( ctable );
		return NULL;
	}
	
	// read file in 2 parts
	ctable->data[0] = ml;
	stream->read( ctable->data + 1, 1, 1 );
	if ( ( ctable->data[0] == 0xFF ) && ( ctable->data[1] == 0xFF ) ) { // special case (easy)
		ctable->max_clen = 8;
		ctable->n_symbols = 256;
		ctable->total_size = 2;
		ctable->len_ptr = NULL;
		ctable->val_ptr = NULL;
	} else { // standard case (read table in 2 parts)
		if ( stream->read( ctable->data + 2, 1, ml - 1 ) != ml - 1 ) {
			free( ctable->data );
			free( ctable );
			return NULL;
		}
		for ( nn = 0, t = 0; t < ml; nn += ctable->data[++t] );
		if ( nn > 256 ) {
			free( ctable->data );
			free( ctable );
			return NULL;
		}
		ctable->data = (unsigned char*) realloc( ctable->data, ml + nn + 1 );
		if ( ctable->data == NULL ) return NULL;
		if ( stream->read( ctable->data + 1 + ml, 1, nn ) != nn ) {
			free( ctable->data );
			free( ctable );
			return NULL;
		}
		ctable->max_clen = ml;
		ctable->n_symbols = nn;
		ctable->total_size = 1 + ml + nn;
		ctable->len_ptr = ctable->data + 1;
		ctable->val_ptr = ctable->data + 1 + ml;
	}
	
	
	return ctable;
}


/* -----------------------------------------------
	free memory from ctable
	----------------------------------------------- */	
void cleanup_ctable( canonical_table* ctable ) {
	if ( ctable != NULL ) {
		if ( ctable->counts != NULL ) free( ctable->counts );
		if ( ctable->data != NULL ) free( ctable->data );
		free( ctable );
	}
}


/* -----------------------------------------------
	quicksort huffman nodes by counts (helper)
	----------------------------------------------- */	
void quicksort_nodes( huffman_node** lbound, huffman_node** rbound ) {
	huffman_node** l;
	huffman_node** r;
	huffman_node* s;
	unsigned int p;
	
	
	// finished if this happens
	if ( lbound >= rbound ) return;
	
	// --- step 1: divide... ---
	l = lbound;
	r = rbound - 1;
	p = (*rbound)->count;
	
	// search for new pivot
	while ( true ) {
		// from left: search element smaller than pivot
		for ( ; ( (*l)->count >= p ) && ( l < rbound ); l++ );
		// from right: search element bigger than pivot
		for ( ; ( (*r)->count <= p ) && ( r > lbound ); r-- );
		// swap if l still on the left, r on the right
		if ( l < r ) {
			s = (*l); (*l) = (*r); (*r) = s;
		} else break; // ...break otherwise
	}
	
	// swap pivot data
	if ( (*l)->count < p ) {
		s = (*l); (*l) = (*rbound); (*rbound) = s;
	}
	
	// --- step 2: ... and conquer! ---
	quicksort_nodes( lbound, l - 1 );
	quicksort_nodes( l + 1, rbound );
	
	
	return;
}


/* -----------------------------------------------
	count code lengths of end nodes (helper)
	----------------------------------------------- */
static void count_lengths( huffman_node* node, unsigned char* clen ) {
	if ( ( node->next[0] == NULL ) && ( node->next[1] == NULL ) ) (*clen)++;
	else {
		clen++;
		if ( node->next[0] != NULL ) count_lengths( node->next[0], clen );
		if ( node->next[1] != NULL ) count_lengths( node->next[1], clen );
	}
}


/* -----------------------------------------------
	free all nodes in tree (helper)
	----------------------------------------------- */	
static void cleanup_tree( huffman_node* node ) {
	if ( node->next[0] != NULL ) cleanup_tree( node->next[0] );
	if ( node->next[1] != NULL ) cleanup_tree( node->next[1] );
	free( node );
}


/* -----------------------------------------------
	convert canonical to hcode structs
	----------------------------------------------- */
	
huffman_code** convert_ctable2hcode( canonical_table* ctable ) {
	huffman_code** table;
	huffman_code* tcode;
	unsigned char* value;
	unsigned char* nlen;
	unsigned int code;
	int ncodes;
	int len;
	int nl;
	
	
	// safety check
	if ( ctable == NULL ) return NULL;
	
	// build and preset huffman code table
	table = ( huffman_code** ) calloc( 256, sizeof( huffman_code* ) );
	if ( table == NULL ) return NULL;
	for ( int i = 0; i < 256; i++ ) table[i] = NULL;
	
	// setup tables
	if ( ctable->len_ptr == NULL ) { // incompressible case table setup
		// build codes - 1 to 1 symbol representation
		for ( code = 0; code < 256; code++ ) {
			tcode = ( huffman_code* ) calloc( 1, sizeof( huffman_code ) );
			if ( tcode == NULL ) { free( table ); return NULL; }
			tcode->len = 8;
			tcode->code = code;
			table[ code ] = tcode;
		}
	} else if ( ctable->n_symbols == 1 ) { // only one symbol in the table
		tcode = ( huffman_code* ) calloc( 1, sizeof( huffman_code ) );
		if ( tcode == NULL ) { free( table ); return NULL; }
		tcode->len = 0;
		tcode->code = 0;
		table[ *(ctable->val_ptr) ] = tcode;
	} else { // standard table setup
		// set up canonical table pointers
		nlen = ctable->len_ptr;
		value = ctable->val_ptr;
		ncodes = ctable->n_symbols;
		
		// build actual codes
		for ( len = 1, code = 0; ncodes > 0; len++ ) {
			for ( nl = (*nlen++); nl > 0; nl--, ncodes-- ) {
				tcode = ( huffman_code* ) calloc( 1, sizeof( huffman_code ) );
				if ( tcode == NULL ) { free( table ); return NULL; }
				tcode->len = len;
				tcode->code = code++;
				table[ (*value++) ] = tcode;
			}
			code <<= 1;
		}
	}
	
	
	return table;
}


/* -----------------------------------------------
	convert canonical to hconv recursive helper
	----------------------------------------------- */
	
huffman_conv_set* convert_ctable2hconv_rec( int* len_ptr, int** val_ptr, int mlen ) {
	huffman_conv_set* conv_set;
	huffman_conv** table;
	huffman_conv* tconv;
	unsigned int c0, c1;
	unsigned int r;
	int hwidth_l;
	int len;	
	int t;
	
	
	// find local hwidth and range
	if ( mlen <= MAX_HWIDTH ) hwidth_l = mlen;
	else for ( r = 1 << MAX_HWIDTH, t = MAX_HWIDTH, hwidth_l = 0;
		( r > 0 ) && ( hwidth_l < MAX_HWIDTH );
		r -= len_ptr[ hwidth_l++ ] * ( 1 << (--t) ) );
	r = 1 << hwidth_l;
	
	// found, allocate memory
	conv_set = ( huffman_conv_set* ) calloc( 1, sizeof( huffman_conv_set ) );
	table = ( huffman_conv** ) calloc( r, sizeof( huffman_conv* ) );
	if ( ( conv_set == NULL ) || ( table == NULL ) ) {
		if ( conv_set != NULL ) free( conv_set );
		if ( table != NULL ) free( table );
		return NULL;
	}
	
	// basic set up of conv_set struct
	for ( c0 = 0; c0 < r; c0++ ) table[ c0 ] = NULL;
	conv_set->hwidth = hwidth_l;
	conv_set->h = table;
	
	// actual generation of conversion table
	for ( len = 1, c0 = c1 = 0; c0 < r; ) {
		// find length
		for ( ; len_ptr[ len ] == 0; len++ );
		// generate new entry
		tconv = ( huffman_conv* ) calloc( 1, sizeof( huffman_conv ) );
		if ( tconv == NULL ) { free( conv_set ); return NULL; }
		
		if ( len <= hwidth_l ) { // routine for regular conversion table
			// set conversion attributes
			tconv->len = len;
			tconv->val = (**val_ptr);
			tconv->ext = NULL;
			// insert it into the table
			for ( c1 += 1 << ( hwidth_l - len ); c0 < c1; c0++ )
				table[ c0 ] = tconv;
			len_ptr[ len ]--;
			(*val_ptr)++;
		} else { // routine for ext conversion table
			// set conversion attributes
			tconv->len = hwidth_l;
			tconv->val = -1;
			tconv->ext = convert_ctable2hconv_rec( len_ptr + hwidth_l, val_ptr, mlen - hwidth_l );
			// insert it into the table
			for ( c1++; c0 < c1; c0++ )
				table[ c0 ] = tconv;
		}
	}
	
	
	return conv_set;
}
	
/* -----------------------------------------------
	convert canonical to hconv structs
	----------------------------------------------- */
	
huffman_conv_set* convert_ctable2hconv( canonical_table* ctable ) {
	huffman_conv_set* conv_set;
	int* llen;
	int* lval;
	int* vptr;
	
	
	// special case: only one symbol in table
	if ( ctable->n_symbols == 1 ) {
		conv_set = ( huffman_conv_set* ) calloc( 1, sizeof( huffman_conv_set ) );
		if ( conv_set == NULL ) return NULL;
		conv_set->h = ( huffman_conv** ) calloc( 1, sizeof( huffman_conv* ) );
		if ( conv_set->h == NULL ) { free( conv_set ); return NULL; }
		conv_set->h[0] = ( huffman_conv* ) calloc( 1, sizeof( huffman_conv ) );
		if ( conv_set->h[0] == NULL ) {	free( conv_set->h[0] );
			free( conv_set->h ); free( conv_set ); return NULL; }
		conv_set->h[0]->len = 0;
		conv_set->h[0]->val = (int) ctable->val_ptr[0]; 
		conv_set->h[0]->ext = NULL;
		conv_set->hwidth = 0;
		return conv_set;
	}
	
	// make a local copy of values and lengths array
	llen = (int*) calloc( ctable->max_clen + 1, sizeof( int ) );
	lval = (int*) calloc( ctable->n_symbols, sizeof( int ) );
	if ( ( llen == NULL ) || ( lval == NULL ) ) {
		if ( llen != NULL ) free( llen );
		if ( lval != NULL ) free( lval );
		return NULL; 
	}
	if ( ctable->len_ptr != NULL ) { // standard case
		for ( int i = 0; i < ctable->max_clen; i++ ) llen[ i + 1 ] = (int) ctable->len_ptr[ i ];
		for ( int i = 0; i < ctable->n_symbols; i++ ) lval[ i ] = (int) ctable->val_ptr[ i ];
	} else { // incompressible case
		llen[ 8 ] = 256;
		for ( int i = 0; i < 256; i++ ) lval[ i ] = i;
	}
	vptr = lval;
	
	// all the heavy lifing is done in a recursive function
	conv_set = convert_ctable2hconv_rec( llen, &vptr, ctable->max_clen );
	
	// free memory
	free( llen );
	free( lval );
	
	
	return conv_set;
}


/* -----------------------------------------------
	cleanup hcode table
	----------------------------------------------- */
void cleanup_hcode( huffman_code** table ) {
	for ( int i = 0; i < 256; i++ )
		if ( table[i] != NULL ) free( table[i] );
	free( table );
}


/* -----------------------------------------------
	cleanup hconv table
	----------------------------------------------- */
void cleanup_hconv( huffman_conv_set* conv_set ) {
	huffman_conv* tconv;
	huffman_conv** table = conv_set->h;
	int r = 1 << conv_set->hwidth;
	int c = 1;
	
	for ( tconv = table[0]; c < r; c++ ) {
		if ( tconv != table[c] ) {
			if ( tconv->ext != NULL ) cleanup_hconv( (huffman_conv_set*) tconv->ext );
			free( tconv );
			tconv = table[c];
		}			
	}
	free( tconv );
	free( table );
	free( conv_set );
}


/* -----------------------------------------------
	constructor for huffman reader class
	----------------------------------------------- */
	
huffman_reader::huffman_reader( unsigned char* data, int size )
{
	// WARNING: NO error checks in this class!
	// init bitreader
	bit_reader = new abitreader( data, size );
	// fill the bit buffer for the first time
	bit_buffer = bit_reader->read( BIT_BUFFER_SIZE );
}


/* -----------------------------------------------
	destructor for huffman reader class
	----------------------------------------------- */
	
huffman_reader::~huffman_reader( void )
{
	// close bitreader
	delete( bit_reader );
}


/* -----------------------------------------------
	decode one symbol from bitstream
	----------------------------------------------- */
	
int huffman_reader::decode_symbol( huffman_conv_set* conv_set )
{
	huffman_conv* conv;
	
	
	do {
		conv = conv_set->h[ bit_buffer >> ( BIT_BUFFER_SIZE - conv_set->hwidth ) ];
		conv_set = (huffman_conv_set*) conv->ext;
		advance_bitstream( conv->len );
	} while ( conv_set != NULL );
	
	
	return conv->val;
}


/* -----------------------------------------------
	bit reader function for n bit
	----------------------------------------------- */
	
unsigned int huffman_reader::read_bits( int n )
{
	unsigned int bits;
	
	// read bits
	bits = ( n <= BIT_BUFFER_SIZE ) ?
		bit_buffer >> (BIT_BUFFER_SIZE-n) :
		( bit_buffer << (n-BIT_BUFFER_SIZE) ) | bit_reader->read( (n-BIT_BUFFER_SIZE) );
	// refill the buffer
	advance_bitstream( n );
	
	return bits;
}


/* -----------------------------------------------
	bit reader function for one bit
	----------------------------------------------- */
	
unsigned char huffman_reader::read_bit( void )
{
	unsigned char bit;
	
	// read one bit
	bit = bit_buffer >> (BIT_BUFFER_SIZE-1);
	// refill the buffer
	advance_bitstream_1();
	
	return bit;
}


/* -----------------------------------------------
	return current byte position
	----------------------------------------------- */
	
int  huffman_reader::getpos( void )
{
	// ugly, but it works!
	return bit_reader->getpos() - ( BIT_BUFFER_SIZE + ( 8 - bit_reader->getbitp() ) ) / 8;
}


/* -----------------------------------------------
	refill bit buffer utility function
	----------------------------------------------- */
	
inline void huffman_reader::advance_bitstream( int n )
{
	// refill the buffer
	bit_buffer = ( n >= BIT_BUFFER_SIZE ) ?
		bit_reader->read( BIT_BUFFER_SIZE ) :
		( ( bit_buffer << n ) | bit_reader->read( n ) ) & ( ( 1 << BIT_BUFFER_SIZE ) - 1 );
}


/* -----------------------------------------------
	refill bit buffer utility function (1 bit)
	----------------------------------------------- */
	
inline void huffman_reader::advance_bitstream_1( void )
{
	// refill the buffer
	bit_buffer = ( ( bit_buffer << 1 ) | bit_reader->read_bit() ) &
		( ( 1 << BIT_BUFFER_SIZE ) - 1 );
}


/* -----------------------------------------------
	constructor for huffman writer class
	----------------------------------------------- */
	
huffman_writer::huffman_writer( int adds )
{
	// WARNING: NO error checks in this class!
	// init bitwriter - recommended value: 5MB
	if ( adds == 0 ) adds = 5 * 1024 * 1024;
	bit_writer = new abitwriter( adds );
}


/* -----------------------------------------------
	destructor for huffman writer class
	----------------------------------------------- */
	
huffman_writer::~huffman_writer( void )
{
	// close bitwriter
	delete( bit_writer );
}


/* -----------------------------------------------
	encode one symbol to bitstream
	----------------------------------------------- */
	
void huffman_writer::encode_symbol( huffman_code** hcodes, int symbol )
{
	huffman_code* hcode;	
	
	// find correct code, encode to bitstream
	hcode = hcodes[symbol];
	write_bits( hcode->code, hcode->len );
}


/* -----------------------------------------------
	bit writer function for n bit
	----------------------------------------------- */
	
void huffman_writer::write_bits( unsigned int bits, int n )
{
	// write bits
	bit_writer->write( bits, n );
}


/* -----------------------------------------------
	bit writer function for 1 bit
	----------------------------------------------- */
	
void huffman_writer::write_bit( unsigned char bit )
{
	// write bit
	bit_writer->write_bit( bit );
}


/* -----------------------------------------------
	return data pointer (when finished)
	----------------------------------------------- */
	
unsigned char* huffman_writer::getptr( void )
{
	return bit_writer->getptr();
}


/* -----------------------------------------------
	return current byte position
	----------------------------------------------- */
	
int huffman_writer::getpos( void )
{
	return bit_writer->getpos();
}
