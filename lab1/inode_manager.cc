#include "inode_manager.h"
/**
 * 30/09/2016
 * editor name: Zhitong Yan
 * id: 5140219099
 */

// disk layer -----------------------------------------
// the layout of the disk should be like this:
// superBLK-freeBLKBitmap-inodeTable-data

disk::disk()
{
  bzero(blocks, sizeof(blocks));
}

void
disk::read_block(blockid_t id, char *buf)
{
  /*
   *your lab1 code goes here.
   *if id is smaller than 0 or larger than BLOCK_NUM 
   *or buf is null, just return.
   *put the content of target block into buf.
   *hint: use memcpy
  */
  	if ((id < 0) || (id > BLOCK_NUM) || (buf == NULL)) {
    		return;
  	}
  	memcpy(buf, blocks[id], BLOCK_SIZE);
  	return;
}

void
disk::write_block(blockid_t id, const char *buf)
{
  /*
   *your lab1 code goes here.
   *hint: just like read_block
  */
  	if ((id < 0) || (id > BLOCK_NUM) || (buf == NULL)) {
    		return;
  	}
  	memcpy(blocks[id], buf, BLOCK_SIZE);
  	return;
}

// block layer -----------------------------------------

// Allocate a free disk block.
// find a free BLK on the disk, and return its blockid
blockid_t
block_manager::alloc_block()
{
  /*
   * your lab1 code goes here.
   * note: you should mark the corresponding bit in block bitmap when alloc.
   * you need to think about which block you can start to be allocated.

   *hint: use macro IBLOCK and BBLOCK.
          use bit operation.
          remind yourself of the layout of disk.
   */
   	
	char block_buf[BLOCK_SIZE];
	blockid_t bid = IBLOCK(sb.ninodes, sb.nblocks);

  	// find the first free BLK number of the free_number_list,
	// and remove from list, and return this number.
	for (uint32_t i = 0; i < BLOCK_NUM; i++) {
		read_block(BBLOCK(bid), block_buf);
		uint32_t index = bid % BLOCK_SIZE;
		uint32_t _index = (index / sizeof(uint32_t));
		uint32_t *mask = &((uint32_t *)block_buf)[_index];
		
		if (!((*mask) & (1 << index))) {
			*mask = *mask | (1 << index);
			write_block(BBLOCK(bid), block_buf);
			bid++;
			return bid-1;	// ??? return bid; 
		} else {
			bid++;
			if (bid >= BLOCK_NUM) {
				bid = bid%BLOCK_NUM + IBLOCK(INODE_NUM,BLOCK_NUM);
			}
		}
	}
  	return 0;
}

void
block_manager::free_block(uint32_t id)
{
  /* 
   * your lab1 code goes here.
   * note: you should unmark the corresponding bit in the block bitmap when free.
   */
	char block_buf[BLOCK_SIZE];
	char mask;
	
	if (id < IBLOCK(INODE_NUM, BLOCK_NUM)-1)
		exit(0);
	read_block(BBLOCK(id), block_buf);
	uint32_t index = id % BPB;
	mask = 1 << (index % 8);
	block_buf[index/8] = block_buf[index/8] & (~mask);
	write_block(BBLOCK(id), block_buf);
	return;  
}

// The layout of disk should be like this:
// |<-sb->|<-free block bitmap->|<-inode table->|<-data->|
block_manager::block_manager()
{
  d = new disk();

  // format the disk
  sb.size = BLOCK_SIZE * BLOCK_NUM;
  sb.nblocks = BLOCK_NUM;
  sb.ninodes = INODE_NUM;

}

void
block_manager::read_block(uint32_t id, char *buf)
{
  d->read_block(id, buf);
}

void
block_manager::write_block(uint32_t id, const char *buf)
{
  d->write_block(id, buf);
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
