#include    <stdint.h>
#include    <stdio.h>
#include    <time.h>

#include        "hardware/vreg.h"
#include        "pico/stdlib.h"
#include	"pico/multicore.h"
#define         LED_PIN         PICO_DEFAULT_LED_PIN


#define MULTITHREAD         1  
#define HAS_FLOAT           1
#define TOTAL_DATA_SIZE     2 * 1000

#define ID_LIST             (1 << 0)
#define ID_MATRIX           (1 << 1)
#define ID_STATE            (1 << 2)
#define ALL_ALGORITHMS_MASK (ID_LIST | ID_MATRIX | ID_STATE)
#define NUM_ALGORITHMS      3

typedef int16_t MATDAT;
typedef int32_t MATRES;
typedef int32_t CORE_TICKS;

typedef struct list_data_s
{
    int16_t data16;
    int16_t idx;
} list_data;

typedef struct list_head_s
{
    struct list_head_s *next;
    struct list_data_s *info;
} list_head;

typedef struct MAT_PARAMS_S
{
    int     N;
    MATDAT *A;
    MATDAT *B;
    MATRES *C;
} mat_params;

typedef struct CORE_PORTABLE_S
{
    uint8_t portable_id;
} core_portable;

typedef struct RESULTS_S
{
    /* inputs */
    int16_t              seed1;       /* Initializing seed */
    int16_t              seed2;       /* Initializing seed */
    int16_t              seed3;       /* Initializing seed */
    void *              memblock[4]; /* Pointer to safe memory location */
    uint32_t              size;        /* Size of the data */
    uint32_t              iterations;  /* Number of iterations to execute */
    uint32_t              execs;       /* Bitmask of operations to execute */
    struct list_head_s *list;
    mat_params          mat;
    /* outputs */
    uint16_t crc;
    uint16_t crclist;
    uint16_t crcmatrix;
    uint16_t crcstate;
    int16_t err;
    /* multithread specific */
    core_portable port;
} core_results;

typedef enum CORE_STATE
{
    CORE_START = 0,
    CORE_INVALID,
    CORE_S1,
    CORE_S2,
    CORE_INT,
    CORE_FLOAT,
    CORE_EXPONENT,
    CORE_SCIENTIFIC,
    NUM_CORE_STATES
} core_state_e;

typedef int32_t (*list_cmp)(list_data *a, list_data *b, core_results *res);

#define matrix_test_next(x)      (x + 1)
#define matrix_clip(x, y)        ((y) ? (x)&0x0ff : (x)&0x0ffff)
#define matrix_big(x)            (0xf000 | (x))
#define bit_extract(x, from, to) (((x) >> (from)) & (~(0xffffffff << (to))))
#define get_seed(x) (int16_t) get_seed_32(x)
#define align_mem(x) (void *)(4 + (((uint32_t)(x)-1) & ~3))

int16_t calc_func(int16_t *pdata, core_results *res);

uint16_t
crcu8(uint8_t data, uint16_t crc)
{
    uint8_t i = 0, x16 = 0, carry = 0;

    for (i = 0; i < 8; i++)
    {
        x16 = (uint8_t)((data & 1) ^ ((uint8_t)crc & 1));
        data >>= 1;

        if (x16 == 1)
        {
            crc ^= 0x4002;
            carry = 1;
        }
        else
            carry = 0;
        crc >>= 1;
        if (carry)
            crc |= 0x8000;
        else
            crc &= 0x7fff;
    }
    return crc;
}
uint16_t
crcu16(uint16_t newval, uint16_t crc)
{
    crc = crcu8((uint8_t)(newval), crc);
    crc = crcu8((uint8_t)((newval) >> 8), crc);
    return crc;
}
uint16_t
crc16(int16_t newval, uint16_t crc)
{
    return crcu16((uint16_t)newval, crc);
}
uint16_t
crcu32(uint32_t newval, uint16_t crc)
{
    crc = crc16((int16_t)newval, crc);
    crc = crc16((int16_t)(newval >> 16), crc);
    return crc;
}

uint32_t core_init_matrix(uint32_t blksize, void *memblk, int32_t seed, mat_params *p)
{
    uint32_t  N = 0;
    MATDAT *A;
    MATDAT *B;
    int32_t  order = 1;
    MATDAT  val;
    uint32_t  i = 0, j = 0;
    if (seed == 0)
        seed = 1;
    while (j < blksize)
    {
        i++;
        j = i * i * 2 * 4;
    }
    N = i - 1;
    A = (MATDAT *)align_mem(memblk);
    B = A + N * N;

    for (i = 0; i < N; i++)
    {
        for (j = 0; j < N; j++)
        {
            seed         = ((order * seed) % 65536);
            val          = (seed + order);
            val          = matrix_clip(val, 0);
            B[i * N + j] = val;
            val          = (val + order);
            val          = matrix_clip(val, 1);
            A[i * N + j] = val;
            order++;
        }
    }

    p->A = A;
    p->B = B;
    p->C = (MATRES *)align_mem(B + N * N);
    p->N = N;

    return N;
}

int16_t matrix_sum(uint32_t N, MATRES *C, MATDAT clipval)
{
    MATRES tmp = 0, prev = 0, cur = 0;
    int16_t ret = 0;
    uint32_t i, j;
    for (i = 0; i < N; i++)
    {
        for (j = 0; j < N; j++)
        {
            cur = C[i * N + j];
            tmp += cur;
            if (tmp > clipval)
            {
                ret += 10;
                tmp = 0;
            }
            else
            {
                ret += (cur > prev) ? 1 : 0;
            }
            prev = cur;
        }
    }
    return ret;
}

void matrix_mul_const(uint32_t N, MATRES *C, MATDAT *A, MATDAT val)
{
    uint32_t i, j;
    for (i = 0; i < N; i++)
    {
        for (j = 0; j < N; j++)
        {
            C[i * N + j] = (MATRES)A[i * N + j] * (MATRES)val;
        }
    }
}

void matrix_add_const(uint32_t N, MATDAT *A, MATDAT val)
{
    uint32_t i, j;
    for (i = 0; i < N; i++)
    {
        for (j = 0; j < N; j++)
        {
            A[i * N + j] += val;
        }
    }
}

void matrix_mul_vect(uint32_t N, MATRES *C, MATDAT *A, MATDAT *B)
{
    uint32_t i, j;
    for (i = 0; i < N; i++)
    {
        C[i] = 0;
        for (j = 0; j < N; j++)
        {
            C[i] += (MATRES)A[i * N + j] * (MATRES)B[j];
        }
    }
}

void matrix_mul_matrix(uint32_t N, MATRES *C, MATDAT *A, MATDAT *B)
{
    uint32_t i, j, k;
    for (i = 0; i < N; i++)
    {
        for (j = 0; j < N; j++)
        {
            C[i * N + j] = 0;
            for (k = 0; k < N; k++)
            {
                C[i * N + j] += (MATRES)A[i * N + k] * (MATRES)B[k * N + j];
            }
        }
    }
}

void matrix_mul_matrix_bitextract(uint32_t N, MATRES *C, MATDAT *A, MATDAT *B)
{
    uint32_t i, j, k;
    for (i = 0; i < N; i++)
    {
        for (j = 0; j < N; j++)
        {
            C[i * N + j] = 0;
            for (k = 0; k < N; k++)
            {
                MATRES tmp = (MATRES)A[i * N + k] * (MATRES)B[k * N + j];
                C[i * N + j] += bit_extract(tmp, 2, 4) * bit_extract(tmp, 5, 7);
            }
        }
    }
}

int16_t matrix_test(uint32_t N, MATRES *C, MATDAT *A, MATDAT *B, MATDAT val)
{
    uint16_t crc     = 0;
    MATDAT clipval = matrix_big(val);

    matrix_add_const(N, A, val); /* make sure data changes  */
    matrix_mul_const(N, C, A, val);
    crc = crc16(matrix_sum(N, C, clipval), crc);

    matrix_mul_vect(N, C, A, B);
    crc = crc16(matrix_sum(N, C, clipval), crc);

    matrix_mul_matrix(N, C, A, B);
    crc = crc16(matrix_sum(N, C, clipval), crc);

    matrix_mul_matrix_bitextract(N, C, A, B);
    crc = crc16(matrix_sum(N, C, clipval), crc);

    matrix_add_const(N, A, -val); /* return matrix to initial value */
    return crc;
}

uint16_t core_bench_matrix(mat_params *p, int16_t seed, uint16_t crc)
{
    uint32_t  N   = p->N;
    MATRES *C   = p->C;
    MATDAT *A   = p->A;
    MATDAT *B   = p->B;
    MATDAT  val = (MATDAT)seed;

    crc = crc16(matrix_test(N, C, A, B, val), crc);

    return crc;
}

int32_t cmp_idx(list_data *a, list_data *b, core_results *res)
{
    if (res == NULL)
    {
        a->data16 = (a->data16 & 0xff00) | (0x00ff & (a->data16 >> 8));
        b->data16 = (b->data16 & 0xff00) | (0x00ff & (b->data16 >> 8));
    }
    return a->idx - b->idx;
}

int32_t cmp_complex(list_data *a, list_data *b, core_results *res)
{
    int16_t val1 = calc_func(&(a->data16), res);
    int16_t val2 = calc_func(&(b->data16), res);
    return val1 - val2;
}

void copy_info(list_data *to, list_data *from)
{
    to->data16 = from->data16;
    to->idx    = from->idx;
}

list_head* core_list_insert_new(list_head * insert_point,
                     list_data * info,
                     list_head **memblock,
                     list_data **datablock,
                     list_head * memblock_end,
                     list_data * datablock_end)
{
    list_head *newitem;

    if ((*memblock + 1) >= memblock_end)
        return NULL;
    if ((*datablock + 1) >= datablock_end)
        return NULL;

    newitem = *memblock;
    (*memblock)++;
    newitem->next      = insert_point->next;
    insert_point->next = newitem;

    newitem->info = *datablock;
    (*datablock)++;
    copy_info(newitem->info, info);

    return newitem;
}

list_head* core_list_remove(list_head *item)
{
    list_data *tmp;
    list_head *ret = item->next;
    /* swap data pointers */
    tmp        = item->info;
    item->info = ret->info;
    ret->info  = tmp;
    /* and eliminate item */
    item->next = item->next->next;
    ret->next  = NULL;
    return ret;
}

list_head* core_list_undo_remove(list_head *item_removed, list_head *item_modified)
{
    list_data *tmp;
    /* swap data pointers */
    tmp                 = item_removed->info;
    item_removed->info  = item_modified->info;
    item_modified->info = tmp;
    /* and insert item */
    item_removed->next  = item_modified->next;
    item_modified->next = item_removed;
    return item_removed;
}

list_head* core_list_find(list_head *list, list_data *info)
{
    if (info->idx >= 0)
    {
        while (list && (list->info->idx != info->idx))
            list = list->next;
        return list;
    }
    else
    {
        while (list && ((list->info->data16 & 0xff) != info->data16))
            list = list->next;
        return list;
    }
}

list_head* core_list_reverse(list_head *list)
{
    list_head *next = NULL, *tmp;
    while (list)
    {
        tmp        = list->next;
        list->next = next;
        next       = list;
        list       = tmp;
    }
    return next;
}

list_head* core_list_mergesort(list_head *list, list_cmp cmp, core_results *res)
{
    list_head *p, *q, *e, *tail;
    int32_t     insize, nmerges, psize, qsize, i;

    insize = 1;

    while (1)
    {
        p    = list;
        list = NULL;
        tail = NULL;

        nmerges = 0; /* count number of merges we do in this pass */

        while (p)
        {
            nmerges++; /* there exists a merge to be done */
            /* step `insize' places along from p */
            q     = p;
            psize = 0;
            for (i = 0; i < insize; i++)
            {
                psize++;
                q = q->next;
                if (!q)
                    break;
            }

            /* if q hasn't fallen off end, we have two lists to merge */
            qsize = insize;

            /* now we have two lists; merge them */
            while (psize > 0 || (qsize > 0 && q))
            {

                /* decide whether next element of merge comes from p or q */
                if (psize == 0)
                {
                    /* p is empty; e must come from q. */
                    e = q;
                    q = q->next;
                    qsize--;
                }
                else if (qsize == 0 || !q)
                {
                    /* q is empty; e must come from p. */
                    e = p;
                    p = p->next;
                    psize--;
                }
                else if (cmp(p->info, q->info, res) <= 0)
                {
                    /* First element of p is lower (or same); e must come from
                     * p. */
                    e = p;
                    p = p->next;
                    psize--;
                }
                else
                {
                    /* First element of q is lower; e must come from q. */
                    e = q;
                    q = q->next;
                    qsize--;
                }

                /* add the next element to the merged list */
                if (tail)
                {
                    tail->next = e;
                }
                else
                {
                    list = e;
                }
                tail = e;
            }

            /* now p has stepped `insize' places along, and q has too */
            p = q;
        }

        tail->next = NULL;

        /* If we have done only one merge, we're finished. */
        if (nmerges <= 1) /* allow for nmerges==0, the empty list case */
            return list;

        /* Otherwise repeat, merging lists twice the size */
        insize *= 2;
    }
    return list;
}

uint16_t core_bench_list(core_results *res, int16_t finder_idx)
{
    uint16_t     retval = 0;
    uint16_t     found = 0, missed = 0;
    list_head *list     = res->list;
    int16_t     find_num = res->seed3;
    list_head *this_find;
    list_head *finder, *remover;
    list_data  info;
    int16_t     i;

    info.idx = finder_idx;
    for (i = 0; i < find_num; i++)
    {
        info.data16 = (i & 0xff);
        this_find   = core_list_find(list, &info);
        list        = core_list_reverse(list);
        if (this_find == NULL)
        {
            missed++;
            retval += (list->next->info->data16 >> 8) & 1;
        }
        else
        {
            found++;
            if (this_find->info->data16 & 0x1) /* use found value */
                retval += (this_find->info->data16 >> 9) & 1;
            /* and cache next item at the head of the list (if any) */
            if (this_find->next != NULL)
            {
                finder          = this_find->next;
                this_find->next = finder->next;
                finder->next    = list->next;
                list->next      = finder;
            }
        }
        if (info.idx >= 0)
            info.idx++;
    }
    retval += found * 4 - missed;
    /* sort the list by data content and remove one item*/
    if (finder_idx > 0)
        list = core_list_mergesort(list, cmp_complex, res);
    remover = core_list_remove(list->next);
    /* CRC data content of list from location of index N forward, and then undo
     * remove */
    finder = core_list_find(list, &info);
    if (!finder)
        finder = list->next;
    while (finder)
    {
        retval = crc16(list->info->data16, retval);
        finder = finder->next;
    }
    remover = core_list_undo_remove(remover, list->next);
    /* sort the list by index, in effect returning the list to original state */
    list = core_list_mergesort(list, cmp_idx, NULL);
    /* CRC data content of list */
    finder = list->next;
    while (finder)
    {
        retval = crc16(list->info->data16, retval);
        finder = finder->next;
    }
    return retval;
}

list_head* core_list_init(uint32_t blksize, list_head *memblock, int16_t seed)
{
    /* calculated pointers for the list */
    uint32_t per_item = 16 + sizeof(struct list_data_s);
    uint32_t size     = (blksize / per_item)
                  - 2; /* to accomodate systems with 64b pointers, and make sure
                          same code is executed, set max list elements */
    list_head *memblock_end  = memblock + size;
    list_data *datablock     = (list_data *)(memblock_end);
    list_data *datablock_end = datablock + size;
    /* some useful variables */
    uint32_t     i;
    list_head *finder, *list = memblock;
    list_data  info;

    /* create a fake items for the list head and tail */
    list->next         = NULL;
    list->info         = datablock;
    list->info->idx    = 0x0000;
    list->info->data16 = (int16_t)0x8080;
    memblock++;
    datablock++;
    info.idx    = 0x7fff;
    info.data16 = (int16_t)0xffff;
    core_list_insert_new(
        list, &info, &memblock, &datablock, memblock_end, datablock_end);

    /* then insert size items */
    for (i = 0; i < size; i++)
    {
        uint16_t datpat = ((uint16_t)(seed ^ i) & 0xf);
        uint16_t dat
            = (datpat << 3) | (i & 0x7); /* alternate between algorithms */
        info.data16 = (dat << 8) | dat;  /* fill the data with actual data and
                                            upper bits with rebuild value */
        core_list_insert_new(
            list, &info, &memblock, &datablock, memblock_end, datablock_end);
    }
    /* and now index the list so we know initial seed order of the list */
    finder = list->next;
    i      = 1;
    while (finder->next != NULL)
    {
        if (i < size / 5) /* first 20% of the list in order */
            finder->info->idx = i++;
        else
        {
            uint16_t pat = (uint16_t)(i++ ^ seed); /* get a pseudo random number */
            finder->info->idx = 0x3fff
                                & (((i & 0x07) << 8)
                                   | pat); /* make sure the mixed items end up
                                              after the ones in sequence */
        }
        finder = finder->next;
    }
    list = core_list_mergesort(list, cmp_idx, NULL);
    return list;
}

enum CORE_STATE core_state_transition(uint8_t **instr, uint32_t *transition_count);

uint16_t
core_bench_state(uint32_t blksize,
                 uint8_t *memblock,
                 int16_t seed1,
                 int16_t seed2,
                 int16_t step,
                 uint16_t crc)
{
    uint32_t final_counts[NUM_CORE_STATES];
    uint32_t track_counts[NUM_CORE_STATES];
    uint8_t *p = memblock;
    uint32_t i;

    for (i = 0; i < NUM_CORE_STATES; i++)
    {
        final_counts[i] = track_counts[i] = 0;
    }
    /* run the state machine over the input */
    while (*p != 0)
    {
        enum CORE_STATE fstate = core_state_transition(&p, track_counts);
        final_counts[fstate]++;
    }
    p = memblock;
    while (p < (memblock + blksize))
    { /* insert some corruption */
        if (*p != ',')
            *p ^= (uint8_t)seed1;
        p += step;
    }
    p = memblock;
    /* run the state machine over the input again */
    while (*p != 0)
    {
        enum CORE_STATE fstate = core_state_transition(&p, track_counts);
        final_counts[fstate]++;
    }
    p = memblock;
    while (p < (memblock + blksize))
    { /* undo corruption is seed1 and seed2 are equal */
        if (*p != ',')
            *p ^= (uint8_t)seed2;
        p += step;
    }
    /* end timing */
    for (i = 0; i < NUM_CORE_STATES; i++)
    {
        crc = crcu32(final_counts[i], crc);
        crc = crcu32(track_counts[i], crc);
    }
    return crc;
}

static uint8_t *intpat[4]
    = { (uint8_t *)"5012", (uint8_t *)"1234", (uint8_t *)"-874", (uint8_t *)"+122" };
static uint8_t *floatpat[4] = { (uint8_t *)"35.54400",
                              (uint8_t *)".1234500",
                              (uint8_t *)"-110.700",
                              (uint8_t *)"+0.64400" };
static uint8_t *scipat[4]   = { (uint8_t *)"5.500e+3",
                            (uint8_t *)"-.123e-2",
                            (uint8_t *)"-87e+832",
                            (uint8_t *)"+0.6e-12" };
static uint8_t *errpat[4]   = { (uint8_t *)"T0.3e-1F",
                            (uint8_t *)"-T.T++Tq",
                            (uint8_t *)"1T3.4e4z",
                            (uint8_t *)"34.0e-T^" };

void
core_init_state(uint32_t size, int16_t seed, uint8_t *p)
{
    uint32_t total = 0, next = 0, i;
    uint8_t *buf = 0;
    size--;
    next = 0;
    while ((total + next + 1) < size)
    {
        if (next > 0)
        {
            for (i = 0; i < next; i++)
                *(p + total + i) = buf[i];
            *(p + total + i) = ',';
            total += next + 1;
        }
        seed++;
        switch (seed & 0x7)
        {
            case 0: /* int */
            case 1: /* int */
            case 2: /* int */
                buf  = intpat[(seed >> 3) & 0x3];
                next = 4;
                break;
            case 3: /* float */
            case 4: /* float */
                buf  = floatpat[(seed >> 3) & 0x3];
                next = 8;
                break;
            case 5: /* scientific */
            case 6: /* scientific */
                buf  = scipat[(seed >> 3) & 0x3];
                next = 8;
                break;
            case 7: /* invalid */
                buf  = errpat[(seed >> 3) & 0x3];
                next = 8;
                break;
            default: /* Never happen, just to make some compilers happy */
                break;
        }
    }
    size++;
    while (total < size)
    { /* fill the rest with 0 */
        *(p + total) = 0;
        total++;
    }
}

static uint8_t
ee_isdigit(uint8_t c)
{
    uint8_t retval;
    retval = ((c >= '0') & (c <= '9')) ? 1 : 0;
    return retval;
}

enum CORE_STATE
core_state_transition(uint8_t **instr, uint32_t *transition_count)
{
    uint8_t *         str = *instr;
    uint8_t           NEXT_SYMBOL;
    enum CORE_STATE state = CORE_START;
    for (; *str && state != CORE_INVALID; str++)
    {
        NEXT_SYMBOL = *str;
        if (NEXT_SYMBOL == ',') /* end of this input */
        {
            str++;
            break;
        }
        switch (state)
        {
            case CORE_START:
                if (ee_isdigit(NEXT_SYMBOL))
                {
                    state = CORE_INT;
                }
                else if (NEXT_SYMBOL == '+' || NEXT_SYMBOL == '-')
                {
                    state = CORE_S1;
                }
                else if (NEXT_SYMBOL == '.')
                {
                    state = CORE_FLOAT;
                }
                else
                {
                    state = CORE_INVALID;
                    transition_count[CORE_INVALID]++;
                }
                transition_count[CORE_START]++;
                break;
            case CORE_S1:
                if (ee_isdigit(NEXT_SYMBOL))
                {
                    state = CORE_INT;
                    transition_count[CORE_S1]++;
                }
                else if (NEXT_SYMBOL == '.')
                {
                    state = CORE_FLOAT;
                    transition_count[CORE_S1]++;
                }
                else
                {
                    state = CORE_INVALID;
                    transition_count[CORE_S1]++;
                }
                break;
            case CORE_INT:
                if (NEXT_SYMBOL == '.')
                {
                    state = CORE_FLOAT;
                    transition_count[CORE_INT]++;
                }
                else if (!ee_isdigit(NEXT_SYMBOL))
                {
                    state = CORE_INVALID;
                    transition_count[CORE_INT]++;
                }
                break;
            case CORE_FLOAT:
                if (NEXT_SYMBOL == 'E' || NEXT_SYMBOL == 'e')
                {
                    state = CORE_S2;
                    transition_count[CORE_FLOAT]++;
                }
                else if (!ee_isdigit(NEXT_SYMBOL))
                {
                    state = CORE_INVALID;
                    transition_count[CORE_FLOAT]++;
                }
                break;
            case CORE_S2:
                if (NEXT_SYMBOL == '+' || NEXT_SYMBOL == '-')
                {
                    state = CORE_EXPONENT;
                    transition_count[CORE_S2]++;
                }
                else
                {
                    state = CORE_INVALID;
                    transition_count[CORE_S2]++;
                }
                break;
            case CORE_EXPONENT:
                if (ee_isdigit(NEXT_SYMBOL))
                {
                    state = CORE_SCIENTIFIC;
                    transition_count[CORE_EXPONENT]++;
                }
                else
                {
                    state = CORE_INVALID;
                    transition_count[CORE_EXPONENT]++;
                }
                break;
            case CORE_SCIENTIFIC:
                if (!ee_isdigit(NEXT_SYMBOL))
                {
                    state = CORE_INVALID;
                    transition_count[CORE_INVALID]++;
                }
                break;
            default:
                break;
        }
    }
    *instr = str;
    return state;
}

volatile int32_t seed1_volatile;
volatile int32_t seed2_volatile;
volatile int32_t seed3_volatile;
volatile int32_t seed4_volatile;
volatile int32_t seed5_volatile;
int32_t
get_seed_32(int i)
{
    int32_t retval;
    switch (i)
    {
        case 1:
            retval = seed1_volatile;
            break;
        case 2:
            retval = seed2_volatile;
            break;
        case 3:
            retval = seed3_volatile;
            break;
        case 4:
            retval = seed4_volatile;
            break;
        case 5:
            retval = seed5_volatile;
            break;
        default:
            retval = 0;
            break;
    }
    return retval;
}

uint8_t
check_data_types()
{
    uint8_t retval = 0;
    if (sizeof(uint8_t) != 1)
    {
        printf("ERROR: uint8_t is not an 8b datatype!\n");
        retval++;
    }
    if (sizeof(uint16_t) != 2)
    {
        printf("ERROR: uint16_t is not a 16b datatype!\n");
        retval++;
    }
    if (sizeof(int16_t) != 2)
    {
        printf("ERROR: int16_t is not a 16b datatype!\n");
        retval++;
    }
    if (sizeof(int32_t) != 4)
    {
        printf("ERROR: int32_t is not a 32b datatype!\n");
        retval++;
    }
    if (sizeof(uint32_t) != 4)
    {
        printf("ERROR: uint32_t is not a 32b datatype!\n");
        retval++;
    }
    if (sizeof(uint32_t) != sizeof(int*))
    {
        printf(
            "ERROR: ee_ptr_int is not a datatype that holds an int pointer!\n");
        retval++;
    }
    if (retval > 0)
    {
        printf("ERROR: Please modify the datatypes in core_portme.h!\n");
    }
    return retval;
}

int16_t calc_func(int16_t *pdata, core_results *res)
{
    int16_t data = *pdata;
    int16_t retval;
    uint8_t  optype
        = (data >> 7)
          & 1;  /* bit 7 indicates if the function result has been cached */
    if (optype) /* if cached, use cache */
        return (data & 0x007f);
    else
    {                             /* otherwise calculate and cache the result */
        int16_t flag = data & 0x7; /* bits 0-2 is type of function to perform */
        int16_t dtype
            = ((data >> 3)
               & 0xf);       /* bits 3-6 is specific data for the operation */
        dtype |= dtype << 4; /* replicate the lower 4 bits to get an 8b value */
        switch (flag)
        {
            case 0:
                if (dtype < 0x22) /* set min period for bit corruption */
                    dtype = 0x22;
                retval = core_bench_state(res->size,
                                          res->memblock[3],
                                          res->seed1,
                                          res->seed2,
                                          dtype,
                                          res->crc);
                if (res->crcstate == 0)
                    res->crcstate = retval;
                break;
            case 1:
                retval = core_bench_matrix(&(res->mat), dtype, res->crc);
                if (res->crcmatrix == 0)
                    res->crcmatrix = retval;
                break;
            default:
                retval = data;
                break;
        }
        res->crc = crcu16(retval, res->crc);
        retval &= 0x007f;
        *pdata = (data & 0xff00) | 0x0080 | retval; /* cache the result */
        return retval;
    }
}


static uint16_t list_known_crc[]   = { (uint16_t)0xd4b0,
                                   (uint16_t)0x3340,
                                   (uint16_t)0x6a79,
                                   (uint16_t)0xe714,
                                   (uint16_t)0xe3c1 };
static uint16_t matrix_known_crc[] = { (uint16_t)0xbe52,
                                     (uint16_t)0x1199,
                                     (uint16_t)0x5608,
                                     (uint16_t)0x1fd7,
                                     (uint16_t)0x0747 };
static uint16_t state_known_crc[]  = { (uint16_t)0x5e47,
                                    (uint16_t)0x39bf,
                                    (uint16_t)0xe5a4,
                                    (uint16_t)0x8e3a,
                                    (uint16_t)0x8d84 };

void* iterate(void* pres)
{
    uint32_t        i;
    uint16_t        crc;
    core_results *res        = (core_results *)pres;
    uint32_t        iterations = res->iterations;
    res->crc                 = 0;
    res->crclist             = 0;
    res->crcmatrix           = 0;
    res->crcstate            = 0;

    for (i = 0; i < iterations; i++)
    {
        crc      = core_bench_list(res, 1);
        res->crc = crcu16(crc, res->crc);
        crc      = core_bench_list(res, -1);
        res->crc = crcu16(crc, res->crc);
        if (i == 0)
            res->crclist = res->crc;
    }
    return NULL;
}

int fakemain(void)
{
    uint16_t       i, j = 0, num_algorithms = 0;
    int16_t       known_id = -1, total_errors = 0;
    uint16_t       seedcrc = 0;
    float          total_time;
    core_results results[MULTITHREAD];
    uint8_t stack_memblock[TOTAL_DATA_SIZE * MULTITHREAD];
    uint32_t default_num_contexts = 1;
    uint64_t start_time;

    if (sizeof(struct list_head_s) > 128)
    {
        printf("list_head structure too big for comparable data!\n");
        return 0;
    }
    results[0].seed1      = get_seed(1);
    results[0].seed2      = get_seed(2);
    results[0].seed3      = get_seed(3);
    results[0].iterations = get_seed_32(4);
    results[0].execs = get_seed_32(5);
    if (results[0].execs == 0)
    { /* if not supplied, execute all algorithms */
        results[0].execs = ALL_ALGORITHMS_MASK;
    }
    /* put in some default values based on one seed only for easy testing */
    if ((results[0].seed1 == 0) && (results[0].seed2 == 0)
        && (results[0].seed3 == 0))
    { /* performance run */
        results[0].seed1 = 0;
        results[0].seed2 = 0;
        results[0].seed3 = 0x66;
    }
    if ((results[0].seed1 == 1) && (results[0].seed2 == 0)
        && (results[0].seed3 == 0))
    { /* validation run */
        results[0].seed1 = 0x3415;
        results[0].seed2 = 0x3415;
        results[0].seed3 = 0x66;
    }
for (i = 0; i < MULTITHREAD; i++)
{
    results[i].memblock[0] = stack_memblock + i * TOTAL_DATA_SIZE;
    results[i].size        = TOTAL_DATA_SIZE;
    results[i].seed1       = results[0].seed1;
    results[i].seed2       = results[0].seed2;
    results[i].seed3       = results[0].seed3;
    results[i].err         = 0;
    results[i].execs       = results[0].execs;
}
    for (i = 0; i < NUM_ALGORITHMS; i++)
    {
        if ((1 << (uint32_t)i) & results[0].execs)
            num_algorithms++;
    }
    for (i = 0; i < MULTITHREAD; i++)
        results[i].size = results[i].size / num_algorithms;
    /* Assign pointers */
    for (i = 0; i < NUM_ALGORITHMS; i++)
    {
        uint32_t ctx;
        if ((1 << (uint32_t)i) & results[0].execs)
        {
            for (ctx = 0; ctx < MULTITHREAD; ctx++)
                results[ctx].memblock[i + 1]
                    = (char *)(results[ctx].memblock[0]) + results[0].size * j;
            j++;
        }
    }
    /* call inits */
    for (i = 0; i < MULTITHREAD; i++)
    {
        if (results[i].execs & ID_LIST)
        {
            results[i].list = core_list_init(
                results[0].size, results[i].memblock[1], results[i].seed1);
        }
        if (results[i].execs & ID_MATRIX)
        {
            core_init_matrix(results[0].size,
                             results[i].memblock[2],
                             (int32_t)results[i].seed1
                                 | (((int32_t)results[i].seed2) << 16),
                             &(results[i].mat));
        }
        if (results[i].execs & ID_STATE)
        {
            core_init_state(
                results[0].size, results[i].seed1, results[i].memblock[3]);
        }
    }

    /* automatically determine number of iterations if not set */
    if (results[0].iterations == 0)
    {
        double secs_passed = 0;
        uint32_t   divisor;
        results[0].iterations = 1;
        while (secs_passed < (double)1)
        {
            results[0].iterations *= 10;
            start_time = to_ms_since_boot(get_absolute_time());
            iterate(&results[0]);
            secs_passed = (to_ms_since_boot(get_absolute_time()) - start_time) / 1000;
        }
        /* now we know it executes for at least 1 sec, set actual run time at
         * about 10 secs */
        divisor = (uint32_t)secs_passed;
        if (divisor == 0) /* some machines cast float to int as 0 since this
                             conversion is not defined by ANSI, but we know at
                             least one second passed */
            divisor = 1;
        results[0].iterations *= 1 + 10 / divisor;
    }
    /* perform actual benchmark */
    start_time = to_ms_since_boot(get_absolute_time());
    iterate(&results[0]);
    total_time = (to_ms_since_boot(get_absolute_time()) - start_time) / 1000.0;
    /* get a function of the input to report */
    seedcrc = crc16(results[0].seed1, seedcrc);
    seedcrc = crc16(results[0].seed2, seedcrc);
    seedcrc = crc16(results[0].seed3, seedcrc);
    seedcrc = crc16(results[0].size, seedcrc);
    switch (seedcrc)
    {                /* test known output for common seeds */
        case 0x8a02: /* seed1=0, seed2=0, seed3=0x66, size 2000 per algorithm */
            known_id = 0;
            printf("6k performance run parameters for coremark.\n");
            break;
        case 0x7b05: /*  seed1=0x3415, seed2=0x3415, seed3=0x66, size 2000 per
                        algorithm */
            known_id = 1;
            printf("6k validation run parameters for coremark.\n");
            break;
        case 0x4eaf: /* seed1=0x8, seed2=0x8, seed3=0x8, size 400 per algorithm
                      */
            known_id = 2;
            printf("Profile generation run parameters for coremark.\n");
            break;
        case 0xe9f5: /* seed1=0, seed2=0, seed3=0x66, size 666 per algorithm */
            known_id = 3;
            printf("2K performance run parameters for coremark.\n");
            break;
        case 0x18f2: /*  seed1=0x3415, seed2=0x3415, seed3=0x66, size 666 per
                        algorithm */
            known_id = 4;
            printf("2K validation run parameters for coremark.\n");
            break;
        default:
            total_errors = -1;
            break;
    }
    if (known_id >= 0)
    {
        for (i = 0; i < default_num_contexts; i++)
        {
            results[i].err = 0;
            if ((results[i].execs & ID_LIST)
                && (results[i].crclist != list_known_crc[known_id]))
            {
                printf("[%u]ERROR! list crc 0x%04x - should be 0x%04x\n",
                          i,
                          results[i].crclist,
                          list_known_crc[known_id]);
                results[i].err++;
            }
            if ((results[i].execs & ID_MATRIX)
                && (results[i].crcmatrix != matrix_known_crc[known_id]))
            {
                printf("[%u]ERROR! matrix crc 0x%04x - should be 0x%04x\n",
                          i,
                          results[i].crcmatrix,
                          matrix_known_crc[known_id]);
                results[i].err++;
            }
            if ((results[i].execs & ID_STATE)
                && (results[i].crcstate != state_known_crc[known_id]))
            {
                printf("[%u]ERROR! state crc 0x%04x - should be 0x%04x\n",
                          i,
                          results[i].crcstate,
                          state_known_crc[known_id]);
                results[i].err++;
            }
            total_errors += results[i].err;
        }
    }
    total_errors += check_data_types();
    /* and report results */
    printf("CoreMark Size    : %lu\n", (long unsigned)results[0].size);
    printf("Total ticks      : %lu\n", (long unsigned)total_time);
    printf("Total time (secs): %f\n", total_time);
    if (total_time > 0)
        printf("Iterations/Sec   : %f\n",
                  default_num_contexts * results[0].iterations
                      / total_time);
    if (total_time < 10)
    {
        printf(
            "ERROR! Must execute for at least 10 secs for a valid result!\n");
        total_errors++;
    }

    printf("Iterations       : %lu\n",
              (long unsigned)default_num_contexts * results[0].iterations);
    /* output for verification */
    printf("seedcrc          : 0x%04x\n", seedcrc);
    if (results[0].execs & ID_LIST)
        for (i = 0; i < default_num_contexts; i++)
            printf("[%d]crclist       : 0x%04x\n", i, results[i].crclist);
    if (results[0].execs & ID_MATRIX)
        for (i = 0; i < default_num_contexts; i++)
            printf("[%d]crcmatrix     : 0x%04x\n", i, results[i].crcmatrix);
    if (results[0].execs & ID_STATE)
        for (i = 0; i < default_num_contexts; i++)
            printf("[%d]crcstate      : 0x%04x\n", i, results[i].crcstate);
    for (i = 0; i < default_num_contexts; i++)
        printf("[%d]crcfinal      : 0x%04x\n", i, results[i].crc);
    if (total_errors == 0)
    {
        printf(
            "Correct operation validated. See README.md for run and reporting "
            "rules.\n");
#if HAS_FLOAT
        if (known_id == 3)
        {
            printf("CoreMark 1.0 : %f\n",
                      default_num_contexts * results[0].iterations
                          / total_time);
        }
#endif
    }
    if (total_errors > 0)
        printf("Errors detected\n");
    if (total_errors < 0)
        printf(
            "Cannot validate operation for these seed values, please compare "
            "with results on a known platform.\n");

    return 0;            
}

void main1(void)
{
	int i = 1;
	while(1)
	{
		fakemain();
		printf("Core 1 just now. Iteration count %d.\n", i);
		i = i + 1;
	}
	return;
}

int main(void)
{

    sleep_ms(1000);
    vreg_set_voltage(VREG_VOLTAGE_1_30);
    sleep_ms(1000);
    if(set_sys_clock_khz(392000, false))
    {
        sleep_ms(1000);
        gpio_init(PICO_DEFAULT_LED_PIN);
        gpio_set_dir(LED_PIN, GPIO_OUT);
        gpio_put(LED_PIN, 1);
    }

    stdio_init_all();

    multicore_launch_core1(main1);

    int i = 1;
    while(1)
    {
        fakemain();
	printf("Core 0 just now. Iteration count %d.\n", i);
	i = i + 1;
    }

    return 0;
}
