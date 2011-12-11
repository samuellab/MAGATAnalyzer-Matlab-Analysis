classdef TrackPoint
    properties 
        loc = single([0;0]); %location as a 2x1
        ind = int16(0); %frame number
        area = single(0); %contour area
        cov = single ([0;0;0]); %covariance matrix, c11, c12=c21,c22
        locInFile = 0; %location in file
        et = 0; %elapsed time
    end
    
    methods
        tp = fromFile (tp, fid, loadIm, loadContour, camcalinfo)
    end %static methods
    methods
         dx = minus(tp2, tp1) %loc2 - loc1
         d = distance(tp1, tp2)
    end
    
    methods %constructor, must be in classdef file
        function tp = TrackPoint(varargin)
            switch (nargin) 
                 case 0                    
                 case 1
                    if (isa(varargin{1}, 'TrackPoint'))
                        op = varargin{1};
                        tp.loc = op.loc;
                        tp.ind = op.ind;
                        tp.area = op.area;
                        tp.cov = op.cov;
                        tp.locInFile = op.locInFile;
                        tp.et = op.et;
                    else
                        if (all(size(varargin{1}) == [1 2]))
                            tp.loc = varargin{1};
                            %tp = class(tp, 'trackpoint');
                        else
                            disp ('bad input to trackpoint constructor ');
                            disp (varargin);
                            %tp = class([], 'trackpoint');
                        end                        
                    end
                otherwise
                    if (all(size(varargin{1}) == [1 2]))
                        tp.loc = varargin{1};
                        varargin = {varargin{2:end}};
                    else
                        tp.loc = [varargin{1} varargin{2}];
                        varargin = {varargin{3:end}};
                    end
                    flist = {'ind', 'area', 'cov'};
                    while (~isempty(varargin) && ~isempty(flist))
                        tp.(flist{1}) = varargin{1};
                        varargin = {varargin{2:end}};
                        flist = {flist{2:end}};
                    end
                    %tp = class(tp, 'trackpoint');
            end%switch
        end%trackpoint
 
       
        %end
    end %constructor
end
%{
class TrackPoint {

	 public:

             static const bool messagesOn = true;

            /*TrackPoint (double x, double y, int t)
             *
             * creates a track point with location x,y and frame t
             * the ID code is generated automatically and is guaranteed to be
             * different from all other automatically generated IDs
             */
//            TrackPoint (double x, double y, int t);

            /*TrackPoint (double x, double y, int t, int ID)
             *
             * creates a track point with location x,y and frame t
             * the ID code is passed in by the creator and should
             * be different from any other ID code, or track match may break
             */
  //          TrackPoint (double x, double y, int t, int ID);


            /*TrackPoint (double x, double y, double cov[], int t)
             *
             * creates a track point with location x,y and frame t
             * the ID code is generated automatically and is guaranteed to be
             * different from all other automatically generated IDs
             *
             * cov[] is either the covariance matrix, an array of 4 doubles
             * or NULL
             */
            TrackPoint (double x, double y,  double area, const double cov[], int t);

            /*TrackPoint (double x, double y, int t, int ID)
             *
             * creates a track point with location x,y and frame t
             * the ID code is passed in by the creator and should
             * be different from any other ID code, or track match may break
             *
             * cov[] is either the covariance matrix, an array of 4 doubles
             * or NULL
             */
            TrackPoint (double x, double y,  double area, const double cov[], int t,int ID);

            /*TrackPoint (TrackPoint *pt)
             *
             * copies data in *pt to a new track point
             */
            TrackPoint (const TrackPoint *pt);
            /*returns the location as a cvPoint of doubles
             *
             */
            inline CvPoint2D32f getLocation() {
                return cvPoint2D32f(x,y);
            }
            //finds the location relative to the image origin point x0,y0
            //e.g. if x,y = 100, 80 and x0,y0 = 20,10, then returns 80,70
            inline CvPoint2D32f getLocation (double x0, double y0) {
                return cvPoint2D32f(x-x0, y-y0);
            }
            inline CvPoint2D32f getLocation (CvPoint pt) {
                return cvPoint2D32f(x-pt.x, y-pt.y);
            }
             inline CvPoint2D32f getLocation (CvPoint2D32f pt) {
                return cvPoint2D32f(x-pt.x, y-pt.y);
            }
            inline void setCovariance(const double *cov) {
                memcpy(&(this->cov), cov, 4*sizeof(double));
            }

            inline void getCovariance(double *cov) {
                memcpy(cov, &(this->cov), 4*sizeof(double));
            }

            inline void setArea (double area) {
                this->area = area;
            }

            inline double getArea () {
                return area;
            }

            /* changes the point's location
             *
             */
            void setLocation(double x, double y);


            /*getFrame, setFrame
             *
             * gets/sets frame number
             */
            inline int getFrame() {
		return frameNum;
            }

            void setFrame(int t);

            /*getID()
             *
             * returns ID number; no setID function by design
             */
            inline int getID() {
                return idNum;
            }

           

            /* int toDisk (FILE *f)
             * f is a pointer to a BINARY output file
             * writes the trackpoint to disk in the format
             * int frame, float x, float y
             *
             * nonzero return value indicates an error
             *
             * subclasses should override this function to output additional
             * information
             */
            virtual int toDisk(FILE *f);

            /* sizeOnDisk
             *
             * returns the number of bytes that will be written to file
             * when toDisk is called
             */
            virtual int sizeOnDisk();
            /*  static TrackPoint *fromDisk(FILE *f);
             *  reads int frame, float x, float y from disk (assumes f is
             *  open to a track point) and stuffs them into a new trackpoint
             *
             *  note that the ID# is unique but different from the ID# of
             *  the TP when it was saved
             *
             *  subclasses should override this function to match toDisk
             */
            static TrackPoint *fromDisk(FILE *f);

             /*virtual string saveDescription()
             *
             * provides a description of how the string is stored on disk
             */
            virtual std::string saveDescription();

            virtual inline std::string name() {
                return std::string("TrackPoint");
            }


            virtual void draw(IplImage *dst, bool active, int x0 = 0, int y0 = 0, int ptrad = 1);

            virtual void drawConnected(IplImage *dst, bool active, int x0, int y0, TrackPoint *pt);

            /* distance, distSquared
             * distance to another point; distSquared is faster
             *
             */
            inline const double distSquared (const TrackPoint &pt2) {
                return ((pt2.x - x) * (pt2.x - x) + (pt2.y - y) * (pt2.y - y));
            }

            inline const double distance (const TrackPoint &pt2){
                return sqrt (distSquared(pt2));
            }

            inline double distSquared (const TrackPoint *pt2) {
                return ((pt2->x - x) * (pt2->x - x) + (pt2->y - y) * (pt2->y - y));
            }

            inline double distance (const TrackPoint *pt2){
                return sqrt (distSquared(pt2));
            }

            //finds the angle of vertex a(this)c
            inline double vertexAngle (const TrackPoint *a, const TrackPoint *c) {
                return acos(((x - a->x)*(x - c->x) + (y - a->y)*(y - c->y))/(distance(a)*distance(c)));
            }
            inline double vertexAngle (const TrackPoint &a, const TrackPoint &c) {
                return acos(((x - a.x)*(x - c.x) + (y - a.y)*(y - c.y))/(distance(a)*distance(c)));
            }

            YAML::Emitter& toYAML (YAML::Emitter& out);

            virtual inline uchar getTypeCode() {
                return _id_code;
            }
            virtual void setMessageHandler (communicator *mh) {
                this->mh = mh;
            }
	 protected:

            

        private:
            TrackPoint(); // keep anyone from instantiating a track point with no location
            static const uchar _id_code = 0x01;
 };

 #endif
%}


