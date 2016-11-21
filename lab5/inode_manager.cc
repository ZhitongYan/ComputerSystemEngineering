#include "inode_manager.h"
/**
 * 13/11/2016
 * editor name: Zhitong Yan
 * id: 5140219099
 */

/*******************************************************************************************/
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#define mm 10
#define nn  ((1<<mm) - 1) 
#define tt  100
#define kk  (nn - 2 * tt) 

int pp[mm + 1] = { 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1 };
int alpha_to[nn + 1], index_of[nn + 1], gg[nn - kk + 1];

int recoverDisk[BLOCK_NUM][nn];

void generate_gf()
{
	register int i, mask;

	mask = 1;
	alpha_to[mm] = 0;
	for (i = 0; i<mm; i++)
	{
		alpha_to[i] = mask;
		index_of[alpha_to[i]] = i;
		if (pp[i] != 0)
			alpha_to[mm] ^= mask;
		mask <<= 1;
	}
	index_of[alpha_to[mm]] = mm;
	mask >>= 1;
	for (i = mm + 1; i<nn; i++)
	{
		if (alpha_to[i - 1] >= mask)
			alpha_to[i] = alpha_to[mm] ^ ((alpha_to[i - 1] ^ mask) << 1);
		else alpha_to[i] = alpha_to[i - 1] << 1;
		index_of[alpha_to[i]] = i;
	}
	index_of[0] = -1;

	// test for pp[]
	for (i = 0; i <= nn; i++) {
		int flag = 0;
		for (int j = 0; j <= nn; j++) {
			if (alpha_to[j] == i)
				flag++;
		}
		if (flag != 1) {
			printf("error: pp[] wrong or parameter don't match\n");
			break;
		}
	}
}

void gen_poly()
{
	register int i, j;

	gg[0] = 2;    /* primitive element alpha = 2  for GF(2**mm)  */
	gg[1] = 1;    /* g(x) = (X+alpha) initially */
	for (i = 2; i <= nn - kk; i++)
	{
		gg[i] = 1;
		for (j = i - 1; j>0; j--)
			if (gg[j] != 0)  gg[j] = gg[j - 1] ^ alpha_to[(index_of[gg[j]] + i) % nn];
			else gg[j] = gg[j - 1];
		gg[0] = alpha_to[(index_of[gg[0]] + i) % nn];     /* gg[0] can never be zero */
	}
	/* convert gg[] to index form for quicker encoding */
	for (i = 0; i <= nn - kk; i++)  gg[i] = index_of[gg[i]];
}


void encode_rs(int* data, int datasize, int* bb, int bbsize)
{
	register int i, j;
	int feedback;

	for (i = 0; i<nn - kk; i++)   bb[i] = 0;
	for (i = kk - 1; i >= 0; i--)
	{
		feedback = index_of[data[i] ^ bb[nn - kk - 1]];
		if (feedback != -1)
		{
			for (j = nn - kk - 1; j>0; j--)
				if (gg[j] != -1)
					bb[j] = bb[j - 1] ^ alpha_to[(gg[j] + feedback) % nn];
				else
					bb[j] = bb[j - 1];
			bb[0] = alpha_to[(gg[0] + feedback) % nn];
		}
		else
		{
			for (j = nn - kk - 1; j>0; j--)
				bb[j] = bb[j - 1];
			bb[0] = 0;
		};
	};
};


void decode_rs(int* recd, int recdsize)
{
	register int i, j, u, q;
	int elp[nn - kk + 2][nn - kk], d[nn - kk + 2], l[nn - kk + 2], u_lu[nn - kk + 2], s[nn - kk + 1];
	int count = 0, syn_error = 0, root[tt], loc[tt], z[tt + 1], err[nn], reg[tt + 1];

	/* first form the syndromes */
	for (i = 1; i <= nn - kk; i++)
	{
		s[i] = 0;
		for (j = 0; j<nn; j++)
			if (recd[j] != -1)
				s[i] ^= alpha_to[(recd[j] + i*j) % nn];      /* recd[j] in index form */
															 /* convert syndrome from polynomial form to index form  */
		if (s[i] != 0)  syn_error = 1;        /* set flag if non-zero syndrome => error */
		s[i] = index_of[s[i]];
	};

	if (syn_error)       /* if errors, try and correct */
	{
		/* initialise table entries */
		d[0] = 0;           /* index form */
		d[1] = s[1];        /* index form */
		elp[0][0] = 0;      /* index form */
		elp[1][0] = 1;      /* polynomial form */
		for (i = 1; i<nn - kk; i++)
		{
			elp[0][i] = -1;   /* index form */
			elp[1][i] = 0;   /* polynomial form */
		}
		l[0] = 0;
		l[1] = 0;
		u_lu[0] = -1;
		u_lu[1] = 0;
		u = 0;

		do
		{
			u++;
			if (d[u] == -1)
			{
				l[u + 1] = l[u];
				for (i = 0; i <= l[u]; i++)
				{
					elp[u + 1][i] = elp[u][i];
					elp[u][i] = index_of[elp[u][i]];
				}
			}
			else
				/* search for words with greatest u_lu[q] for which d[q]!=0 */
			{
				q = u - 1;
				while ((d[q] == -1) && (q>0)) q--;
				/* have found first non-zero d[q]  */
				if (q>0)
				{
					j = q;
					do
					{
						j--;
						if ((d[j] != -1) && (u_lu[q]<u_lu[j]))
							q = j;
					} while (j>0);
				};

				/* have now found q such that d[u]!=0 and u_lu[q] is maximum */
				/* store degree of new elp polynomial */
				if (l[u]>l[q] + u - q)  l[u + 1] = l[u];
				else  l[u + 1] = l[q] + u - q;

				/* form new elp(x) */
				for (i = 0; i<nn - kk; i++)    elp[u + 1][i] = 0;
				for (i = 0; i <= l[q]; i++)
					if (elp[q][i] != -1)
						elp[u + 1][i + u - q] = alpha_to[(d[u] + nn - d[q] + elp[q][i]) % nn];
				for (i = 0; i <= l[u]; i++)
				{
					elp[u + 1][i] ^= elp[u][i];
					elp[u][i] = index_of[elp[u][i]];  /*convert old elp value to index*/
				}
			}
			u_lu[u + 1] = u - l[u + 1];

			/* form (u+1)th discrepancy */
			if (u<nn - kk)    /* no discrepancy computed on last iteration */
			{
				if (s[u + 1] != -1)
					d[u + 1] = alpha_to[s[u + 1]];
				else
					d[u + 1] = 0;
				for (i = 1; i <= l[u + 1]; i++)
					if ((s[u + 1 - i] != -1) && (elp[u + 1][i] != 0))
						d[u + 1] ^= alpha_to[(s[u + 1 - i] + index_of[elp[u + 1][i]]) % nn];
				d[u + 1] = index_of[d[u + 1]];    /* put d[u+1] into index form */
			}
		} while ((u<nn - kk) && (l[u + 1] <= tt));

		u++;
		if (l[u] <= tt)         /* can correct error */
		{
			/* put elp into index form */
			for (i = 0; i <= l[u]; i++)   elp[u][i] = index_of[elp[u][i]];

			/* find roots of the error location polynomial */
			for (i = 1; i <= l[u]; i++)
				reg[i] = elp[u][i];
			count = 0;
			for (i = 1; i <= nn; i++)
			{
				q = 1;
				for (j = 1; j <= l[u]; j++)
					if (reg[j] != -1)
					{
						reg[j] = (reg[j] + j) % nn;
						q ^= alpha_to[reg[j]];
					};
				if (!q)        /* store root and error location number indices */
				{
					root[count] = i;
					loc[count] = nn - i;
					count++;
				};
			};
			if (count == l[u])    /* no. roots = degree of elp hence <= tt errors */
			{
				/* form polynomial z(x) */
				for (i = 1; i <= l[u]; i++)        /* Z[0] = 1 always - do not need */
				{
					if ((s[i] != -1) && (elp[u][i] != -1))
						z[i] = alpha_to[s[i]] ^ alpha_to[elp[u][i]];
					else if ((s[i] != -1) && (elp[u][i] == -1))
						z[i] = alpha_to[s[i]];
					else if ((s[i] == -1) && (elp[u][i] != -1))
						z[i] = alpha_to[elp[u][i]];
					else
						z[i] = 0;
					for (j = 1; j<i; j++)
						if ((s[j] != -1) && (elp[u][i - j] != -1))
							z[i] ^= alpha_to[(elp[u][i - j] + s[j]) % nn];
					z[i] = index_of[z[i]];         /* put into index form */
				};

				/* evaluate errors at locations given by error location numbers loc[i] */
				for (i = 0; i<nn; i++)
				{
					err[i] = 0;
					if (recd[i] != -1)        /* convert recd[] to polynomial form */
						recd[i] = alpha_to[recd[i]];
					else  recd[i] = 0;
				}
				for (i = 0; i<l[u]; i++)    /* compute numerator of error term first */
				{
					err[loc[i]] = 1;       /* accounts for z[0] */
					for (j = 1; j <= l[u]; j++)
						if (z[j] != -1)
							err[loc[i]] ^= alpha_to[(z[j] + j*root[i]) % nn];
					if (err[loc[i]] != 0)
					{
						err[loc[i]] = index_of[err[loc[i]]];
						q = 0;     /* form denominator of error term */
						for (j = 0; j<l[u]; j++)
							if (j != i)
								q += index_of[1 ^ alpha_to[(loc[j] + root[i]) % nn]];
						q = q % nn;
						err[loc[i]] = alpha_to[(err[loc[i]] - q + nn) % nn];
						recd[loc[i]] ^= err[loc[i]];  /*recd[i] must be in polynomial form */
					}
				}
			}
			else    /* no. roots != degree of elp => >tt errors and cannot solve */
				for (i = 0; i<nn; i++)        /* could return error flag if desired */
					if (recd[i] != -1)        /* convert recd[] to polynomial form */
						recd[i] = alpha_to[recd[i]];
					else  recd[i] = 0;     /* just output received codeword as is */
		}
		else         /* elp has degree has degree >tt hence cannot solve */
			for (i = 0; i<nn; i++)       /* could return error flag if desired */
				if (recd[i] != -1)        /* convert recd[] to polynomial form */
					recd[i] = alpha_to[recd[i]];
				else  recd[i] = 0;     /* just output received codeword as is */
	}
	else       /* no non-zero syndromes => no errors: output received codeword */
		for (i = 0; i<nn; i++)
			if (recd[i] != -1)        /* convert recd[] to polynomial form */
				recd[i] = alpha_to[recd[i]];
			else  recd[i] = 0;
}
/*******************************************************************************************/

unsigned char cpblocks1[BLOCK_NUM][BLOCK_SIZE];

void initcp()
{
	bzero(cpblocks1, sizeof(cpblocks1));
}

void readcp(blockid_t id, char *buf)
{
	if ((id < 0) || (id > BLOCK_NUM) || (buf == NULL)) {
    		return;
  	}
  	memcpy(buf, cpblocks1[id], BLOCK_SIZE);
  	return;
}

void writecp(blockid_t id, const char *buf)
{
  	if ((id < 0) || (id > BLOCK_NUM) || (buf == NULL)) {
    		return;
  	}
  	memcpy(cpblocks1[id], buf, BLOCK_SIZE);
  	return;
}


// disk layer -----------------------------------------
// the layout of the disk should be like this:
// superBLK-freeBLKBitmap-inodeTable-data
// for lab5 comment this part
/*
disk::disk()
{
  bzero(blocks, sizeof(blocks));
}

void
disk::read_block(blockid_t id, char *buf)
{
  	if ((id < 0) || (id > BLOCK_NUM) || (buf == NULL)) {
    		return;
  	}
  	memcpy(buf, blocks[id], BLOCK_SIZE);
  	return;
}

void
disk::write_block(blockid_t id, const char *buf)
{
  	if ((id < 0) || (id > BLOCK_NUM) || (buf == NULL)) {
    		return;
  	}
  	memcpy(blocks[id], buf, BLOCK_SIZE);
  	return;
}
*/

// block layer -----------------------------------------

// this layer was modified because of the bug of lab3...

// Allocate a free disk block.
blockid_t
block_manager::alloc_block()
{
	char buf[BLOCK_SIZE], mask;

	// block number is 1024*32, and BPB=1024*4.
	// 8 bitmap blocks, so i add most 7 times.
	// from the first block! ---- block[0](boot block).
	for (blockid_t i = 0; i < BLOCK_NUM; i += BPB) {
		read_block(BBLOCK(i), buf);
		// buf is the related bitmap block.
		// set the bit
		for (uint32_t offset = 0; (offset < BPB) && (i + offset < BLOCK_NUM); offset++) {
			mask = 1 << (offset % 8);
			uint32_t index = offset / 8;
			if (!(buf[index] & mask)) {
				buf[index] = buf[index] | mask;
				write_block(BBLOCK(i), buf);
				return i + offset;
			}
		}
	}
	exit(0);
}

void
block_manager::free_block(uint32_t id)
{
	char buf[BLOCK_SIZE], mask;
	if (id < 2 + BLOCK_NUM/BPB + INODE_NUM/IPB) {
		// can not free matadata blk.
		exit(0);
	}
	// unset the bit.
	read_block(BBLOCK(id), buf);
  	uint32_t offset = id % BPB;
	mask = 1 << (offset % 8);
	buf[offset/8] &= (~mask);
	write_block(BBLOCK(id), buf);
	return;
}

// The layout of disk should be like this:
// |<-sb->|<-free block bitmap->|<-inode table->|<-data->|
block_manager::block_manager()
{
  d = new disk();
  //initcp();

  	for (int i = 0; i < BLOCK_NUM; i++) {
  		for (int j = 0; j < nn; j++) {
  			recoverDisk[i][j] = 0;
  		}
  	}
  	generate_gf();
  	gen_poly();
  	// pp, alpha_to, index_of, gg
  	// this arrays from here READ ONLY

  // format the disk
  sb.size = BLOCK_SIZE * BLOCK_NUM;
  sb.nblocks = BLOCK_NUM;
  sb.ninodes = INODE_NUM;

	alloc_block();		// boot block
	alloc_block();		// superblock
	for (uint32_t i = 0; i < BLOCK_NUM/BPB; i++) {
		alloc_block();	// bitmap block
	}
	for (uint32_t i = 0; i < INODE_NUM/IPB; i++) {
		alloc_block();	// inode_table block
	}
	// record the super block.
	char buf[BLOCK_SIZE];
	memset(buf, 0, BLOCK_SIZE);
	memcpy(buf, &sb, sizeof(sb));
	write_block(1, buf);
}

// never used...
void readcpy(blockid_t id, char *buf)
{
	int recd[nn];
	for(int i = 0; i < nn; i++) {
		recd[i] = recoverDisk[id][i];
	}
	for (int i = 0; i < nn; i++) {
		recd[i] = index_of[recd[i]];
	}
	decode_rs(recd, nn);

	unsigned char buf2[BLOCK_SIZE];
	for(int i = 0; i < BLOCK_SIZE; i++) {
		buf2[i] = recd[i];
	}
	if ((id < 0) || (id > BLOCK_NUM) || (buf == NULL)) {
    		return;
  	}
  	memcpy(buf, buf2, BLOCK_SIZE);
}

void
block_manager::read_block(uint32_t id, char *buf)
{
	char buf1[BLOCK_SIZE];
	char buf2[BLOCK_SIZE];

	readcp(id, buf1);
	d->read_block(id, buf2);

	if(!strcmp(buf1, buf2))
  		d->read_block(id, buf); // original
  	else
  		readcp(id, buf);
}

void
block_manager::write_block(uint32_t id, const char *buf)
{
	writecp(id, buf);
  	d->write_block(id, buf); // original

  	int data[kk];
  	int bb[nn - kk];
	for (int i = 0; i < kk; i++) {
		data[i] = 257;
	}
	for (int i = 0; i < BLOCK_SIZE; i++) {
		data[i] = buf[i];
	}
	encode_rs(data, kk, bb, nn-kk);

	int recd[nn];
	for (int i = 0; i<nn - kk; i++)  
		recd[i] = bb[i];
	for (int i = 0; i<kk; i++) 
		recd[i + nn - kk] = data[i];

	for (int i = 0; i < nn; i++) {
		recoverDisk[id][i] = recd[i];
	}

	//d->write_block(id, buf);
}

// inode layer -----------------------------------------

inode_manager::inode_manager()
{
  bm = new block_manager();
  uint32_t root_dir = alloc_inode(extent_protocol::T_DIR);
  if (root_dir != 1) {
    printf("\tim: error! alloc first inode %d, should be 1\n", root_dir);
    exit(0);
  }
}

/* Create a new file.
 * Return its inum. */
uint32_t
inode_manager::alloc_inode(uint32_t type)
{
  /* 
   * your lab1 code goes here.
   * note: the normal inode block should begin from the 2nd inode block.
   * the 1st is used for root_dir, see inode_manager::inode_manager().
    
   * if you get some heap memory, do not forget to free it.
   */
	for (uint32_t i = 1; i < INODE_NUM; i++) {
		inode_t *inode_p = get_inode(i);
		if (inode_p == NULL) {
			struct inode _i;
			_i.type = type;
			_i.size = 0;
			_i.atime = 0;
			_i.mtime = 0;
			_i.ctime = 0;
			put_inode(i, &_i);
			return i;
		} else {
			delete[] inode_p;
			continue;
		}
	}
	exit(0);
}

void
inode_manager::free_inode(uint32_t inum)
{
  /* 
   * your lab1 code goes here.
   * note: you need to check if the inode is already a freed one;
   * if not, clear it, and remember to write back to disk.
   * do not forget to free memory if necessary.
   */
	inode_t *inode_p = get_inode(inum);
	if (inode_p != NULL) {
		inode_p->type = 0;
		inode_p->size = 0;
		put_inode(inum, inode_p);
		free(inode_p);
	}
	return;
}


/* Return an inode structure by inum, NULL otherwise.
 * Caller should release the memory. */
struct inode* 
inode_manager::get_inode(uint32_t inum)
{
  struct inode *ino, *ino_disk;
  char buf[BLOCK_SIZE];

  printf("\tim: get_inode %d\n", inum);

  if (inum < 0 || inum >= INODE_NUM) {
    printf("\tim: inum out of range\n");
    return NULL;
  }

  bm->read_block(IBLOCK(inum, bm->sb.nblocks), buf);
  // printf("%s:%d\n", __FILE__, __LINE__);

  ino_disk = (struct inode*)buf + inum%IPB;
  if (ino_disk->type == 0) {
    printf("\tim: inode not exist\n");
    return NULL;
  }

  ino = (struct inode*)malloc(sizeof(struct inode));
  *ino = *ino_disk;

  return ino;
}

void
inode_manager::put_inode(uint32_t inum, struct inode *ino)
{
  char buf[BLOCK_SIZE];
  struct inode *ino_disk;

  printf("\tim: put_inode %d\n", inum);
  if (ino == NULL)
    return;

  bm->read_block(IBLOCK(inum, bm->sb.nblocks), buf);
  ino_disk = (struct inode*)buf + inum%IPB;
  *ino_disk = *ino;
  bm->write_block(IBLOCK(inum, bm->sb.nblocks), buf);
}

#define MIN(a,b) ((a)<(b) ? (a) : (b))

/* Get all the data of a file by inum. 
 * Return alloced data, should be freed by caller. */
void
inode_manager::read_file(uint32_t inum, char **buf_out, int *size)
{
  /*
   * your lab1 code goes here.
   * note: read blocks related to inode number inum,
   * and copy them to buf_out
   */
	inode_t *inode_p = get_inode(inum);
	if (inode_p == NULL) {
		exit(0);
	}
	
	*size = inode_p->size;
	uint32_t bNum = 0;	// the block numbers of this file
	if (!((inode_p->size) % BLOCK_SIZE)) {
		bNum = inode_p->size / BLOCK_SIZE;	
	} else {
		bNum = (inode_p->size / BLOCK_SIZE) + 1;
	}
	*buf_out = (char *)malloc(bNum * BLOCK_SIZE);
	if (bNum <= NDIRECT) {
		for (uint32_t i = 0; i < bNum; i++) {
			bm->read_block(inode_p->blocks[i], *buf_out + i*BLOCK_SIZE);
		}
		delete inode_p;
		return;
	} else {
		// direct blocks
		for (uint32_t i = 0; i < NDIRECT; i++) {
			bm->read_block(inode_p->blocks[i], *buf_out + i*BLOCK_SIZE);
		}
		// indirect block
		//blockid_t indrct_blk[BLOCK_SIZE];
		blockid_t *indrct_blk = (blockid_t *)malloc(BLOCK_SIZE);
		bm->read_block(inode_p->blocks[NDIRECT], (char *)indrct_blk);
		for (uint32_t i = 0; i < bNum-NDIRECT; i++) {
			bm->read_block(indrct_blk[i], *buf_out + (NDIRECT+i)*BLOCK_SIZE);
		}
		delete[] indrct_blk;
		delete inode_p;
		return;
	}
}

/* alloc/free blocks if needed */
void
inode_manager::write_file(uint32_t inum, const char *buf, int size)
{
  /*
   * your lab1 code goes here.
   * note: write buf to blocks of inode inum.
   * you need to consider the situation when the size of buf 
   * is larger or smaller than the size of original inode.
   * you should free some blocks if necessary.
   */
	inode_t *inode_p = get_inode(inum);
	int _indrct_flag = 0;
	if ((inode_p->size/BLOCK_SIZE) > NDIRECT) {
		_indrct_flag = 1;
	}

	// old value
	uint32_t _size = inode_p->size;
	uint32_t _bNum = 0;
	if (!(_size % BLOCK_SIZE)) {
		_bNum = _size / BLOCK_SIZE;
	} else {
		_bNum = _size / BLOCK_SIZE + 1;
	}

	// update
	inode_p->size = MAXFILE*BLOCK_SIZE > (unsigned)size ? (unsigned)size : MAXFILE*BLOCK_SIZE;
	
	uint32_t bNum = 0;
	if (!(inode_p->size % BLOCK_SIZE)) {
		bNum = inode_p->size / BLOCK_SIZE;
	} else {
		bNum = inode_p->size / BLOCK_SIZE + 1;
	}
	int indrct_flag = 0;
	if (bNum > NDIRECT) {
		indrct_flag = 1;
	}
	char indrct_buf[BLOCK_SIZE];
	char w_buf[BLOCK_SIZE];
	uint32_t *indrct_a;
	uint32_t remain_len = inode_p->size;
	uint32_t total_len = 0, w_len = 0, w_blk = 0, w_id = 0;
	
	// free direct blocks
	for (uint32_t i = 0; i < _bNum && i < NDIRECT; i++) {
		bm->free_block(inode_p->blocks[i]);
	}
	// free indirect block
	if (_indrct_flag) {
		uint32_t indrct_blk = inode_p->blocks[NDIRECT];
		char _buf[BLOCK_SIZE];
		bm->read_block(indrct_blk, _buf);
		uint32_t *_indrct_a = (uint32_t *)_buf;
		for (uint32_t i = 0; i < _bNum - NDIRECT; i++) {
			bm->free_block(_indrct_a[i]);
		}
		bm->free_block(indrct_blk);
	}
	
	// init indirect block
	if (indrct_flag) {
		inode_p->blocks[NDIRECT] = bm->alloc_block();
		bm->read_block(inode_p->blocks[NDIRECT], indrct_buf);
		indrct_a = (uint32_t *)indrct_buf;
	}

	// begin write
	for (total_len = 0; w_blk < bNum; total_len += w_len, remain_len -= w_len) {
		w_len = remain_len > BLOCK_SIZE ? BLOCK_SIZE : remain_len;
		// direct
		if (w_blk < NDIRECT) {
			w_id = bm->alloc_block();
			inode_p->blocks[w_blk] = w_id;
		}
		// indirect
		else {
			w_id = bm->alloc_block();
			indrct_a[w_blk - NDIRECT] = w_id;
			// update this indirect blk
			bm->write_block(inode_p->blocks[NDIRECT], indrct_buf);
		}
		memset(w_buf, 0, BLOCK_SIZE);
		memcpy(w_buf, buf, w_len);
		bm->write_block(w_id, w_buf);
		buf += w_len;
		w_blk++;
	}
	// change file matadata
	inode_p->atime = inode_p->ctime = inode_p->mtime = time(NULL);
	put_inode(inum, inode_p);
	free(inode_p);
  	return;
}

void
inode_manager::getattr(uint32_t inum, extent_protocol::attr &a)
{
  /*
   * your lab1 code goes here.
   * note: get the attributes of inode inum.
   * you can refer to "struct attr" in extent_protocol.h
   */
  	inode_t *inode_p = get_inode(inum);
	if (inode_p != NULL) {
    		a.type = inode_p->type;
    		a.atime = inode_p->atime;
    		a.mtime = inode_p->mtime;
    		a.ctime = inode_p->ctime;
     		a.size = inode_p->size;
   	} else {
     		a.type = 0;
     		a.atime = 0;
     		a.mtime = 0;
     		a.ctime = 0;
     		a.size = 0;
   	}
	free(inode_p);
	return;
}

void
inode_manager::remove_file(uint32_t inum)
{
  /*
   * your lab1 code goes here
   * note: you need to consider about both the data block and inode of the file
   * do not forget to free memory if necessary.
   */
	inode_t *inode_p = get_inode(inum);
	if (inode_p == NULL) return;

	uint32_t num = (inode_p->size + BLOCK_SIZE - 1) / BLOCK_SIZE;
	uint32_t min = num > NDIRECT ? NDIRECT : num;
	// remove all blocks related to this file(inode)
	for (uint32_t i = 0; i < min; i++) {
		bm->free_block(inode_p->blocks[i]);	
	}
	
	// indirect blk
	blockid_t id = inode_p->blocks[NDIRECT];
    	if (num > NDIRECT) {
        	char buf[BLOCK_SIZE];
        	bm->read_block(id, buf);	// indirect blk
        	uint32_t *ip = (uint32_t *)buf;
        	for(uint32_t i = 0; i < num - NDIRECT; i++) {
            		bm->free_block(ip[i]);
		}
        	bm->free_block(id);
    	}
    	free(inode_p);
    	free_inode(inum);
}
