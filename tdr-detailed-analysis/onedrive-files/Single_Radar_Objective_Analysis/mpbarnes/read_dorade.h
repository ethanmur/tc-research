/***************************************************/
/* code obtained from Wen-Chau Lee at the National */
/* Center for Atmospheric Research on 8/6/01       */
/***************************************************/

/* modified 2/28/06 by David Dowell:                                                    */
/* -- added compiler option -DLONG32, which reads "int" (32 bit) type instead of "long" */

/***************************************************/
/* CONSTANTS */
/***************************************************/
#define FALSE 0
#define TRUE 1
#define IDENT_LEN 4
#define PROJ_NAME_LEN 20
#define FLIGHT_NUM_LEN 8
#define FAC_NAME_LEN 8
#define RAD_NAME_LEN 8
#define MAX_NUM_PARMS 50 
#define PARM_NAME_LEN 8 
#define PARM_DESC_LEN 40
#define PARM_UNIT_LEN 8
#define THRESHOLD_FLD_LEN 8
#define MAX_BEAMS 1500 
#define MAX_GATES 3000
/***************************************************/
/* STRUCTURES */
/***************************************************/
struct vold_info {

   short ver_num;
   short vol_num;
#ifdef LONG32
   int max_bytes;
#else
   long max_bytes;
#endif
   char proj_name[PROJ_NAME_LEN];
   short year;
   short mon;
   short day;
   short hour;
   short min;
   short sec;
   char flight_num[FLIGHT_NUM_LEN];
   char gen_fac_name[FAC_NAME_LEN];
   short gen_year;
   short gen_mon;
   short gen_day;
   short num_sensors;

};
struct radd_info {
   
   char rad_name[RAD_NAME_LEN];
   float rad_constant;
   float peak_pow;
   float noise_pos;
   float rec_gain;
   float ant_gain;
   float sys_gain;
   float horiz_beam_width;
   float vert_beam_width;
   short rad_type;
   short scan_mode;
   float scan_rate;
   float start_ang;
   float stop_ang;
   short num_param_desc;
   short num_desc;
   short compress_flag;
   short data_reduc_flag;
   float data_reduc_parm1;
   float data_reduc_parm2;
   float radar_lon;
   float radar_lat;
   float radar_alt;
   float unambig_vel;
   float unambig_range;
   short num_freq;
   short num_pulse;
   float freq1;
   float freq2;
   float freq3;
   float freq4;
   float freq5;
   float pulse_per1;
   float pulse_per2;
   float pulse_per3;
   float pulse_per4;
   float pulse_per5;

};

struct cfac_info {
   
   float c_azimuth;
   float c_elevation;
   float c_range_delay;
   float c_rad_lon;
   float c_rad_lat;
   float c_alt_msl;
   float c_alt_agl;
   float c_ew_grspeed;
   float c_ns_grspeed;
   float c_vert_vel;
   float c_head;
   float c_roll;
   float c_pitch;
   float c_drift;
   float c_rotang;
   float c_tiltang;

};

struct parm_info {

   char parm_name[PARM_NAME_LEN];
   char parm_desc[PARM_DESC_LEN];
   char parm_unit[PARM_UNIT_LEN];
   short ipp;
   short trans_freq;
   float rec_bandwidth;
   short pulse_width;
   short polarization;
   short num_samples;
   short parm_type;
   char threshold_fld[THRESHOLD_FLD_LEN];
   float threshold_val;
   float scale_fac;
   float offset_fac;
#ifdef LONG32
   int baddata_flag;
#else
   long baddata_flag;
#endif
    
};

struct celv_info {

#ifdef LONG32
   int total_gates;
#else
   long total_gates;
#endif
   float gate_spacing[MAX_GATES];

};

struct swib_info {

   char rad_name[RAD_NAME_LEN];
#ifdef LONG32
   int sweep_num;
   int num_rays;
#else
   long sweep_num;
   long num_rays;
#endif
   float start_ang;
   float stop_ang;
   float fixed_ang;
#ifdef LONG32
   int filter_flag;
#else
   long filter_flag;
#endif
   
};

struct ryib_info {

#ifdef LONG32
   int sweep_num;
   int julian_day;
#else
   long sweep_num;
   long julian_day;
#endif
   short hour;
   short min;
   short sec;
   short msec;
   float azimuth;
   float elevation;
   float peak_power;
   float scan_rate;
#ifdef LONG32
   int ray_status;
#else
   long ray_status;
#endif

};

struct asib_info {

   float lon;
   float lat;
   float alt_msl;
   float alt_agl;
   float ew_gspeed;
   float ns_gspeed;
   float vert_vel;
   float head;
   float roll;
   float pitch;
   float drift;
   float rot_ang;
   float tilt_ang;
   float ew_horiz_wind;
   float ns_horiz_wind;
   float vert_wind;
   float head_change;
   float pitch_change;

};

struct rdat_info {

   char parm_name[PARM_NAME_LEN];
   float data[MAX_GATES];

};

/* PROTOTYPES */
void sweepread_(char swp_fname[],struct vold_info *vptr,
                struct radd_info *rptr,struct celv_info *cptr,
                struct cfac_info *cfptr,struct parm_info *pptr,
                struct swib_info *sptr,struct ryib_info *ryptr,
                struct asib_info *aptr,struct rdat_info *dptr,
                char fieldname[]);
void read_vold(FILE *fp,struct vold_info *vptr);
void read_radd(FILE *fp,struct radd_info *rptr);
void read_cfac(FILE *fp,struct cfac_info *cptr);
void read_parm(FILE *fp,struct parm_info *pptr);
void read_celv(FILE *fp,struct celv_info *cptr,int desc_length);
void read_swib(FILE *fp,struct swib_info *sptr);
void read_ryib(FILE *fp,struct ryib_info *rptr);
void read_asib(FILE *fp,struct asib_info *aptr);
#ifdef LONG32
void read_rdat(FILE *fp,char fld_name[],int fld_num,
               int desc_len,int *match,short parm_type,
               int beam_count,int total_gates,
               short compression,int baddata_flag,
               float scale_fac,struct rdat_info *dptr);
void read_sh_arr(FILE *fp,int arrsize,int beam_count,
                 int total_gates,short compression,
                 int baddata_flag,float scale_fac,
                 struct rdat_info *dptr);
#else
void read_rdat(FILE *fp,char fld_name[],int fld_num,
               int desc_len,int *match,short parm_type,
               int beam_count,long total_gates,
               short compression,long baddata_flag,
               float scale_fac,struct rdat_info *dptr);
void read_sh_arr(FILE *fp,int arrsize,int beam_count,
                 long total_gates,short compression,
                 long baddata_flag,float scale_fac,
                 struct rdat_info *dptr);
#endif
void read_ch_arr(FILE *fp,int arrsize);
void read_lg_arr(FILE *fp,int arrsize);
void read_fl_arr(FILE *fp,int arrsize);
void get_field(struct parm_info parm[],int num_desc,int *fld_num);
long read_long(FILE *fp);
int read_int(FILE *fp);
void skip_bytes(FILE *fp,int numskip);
