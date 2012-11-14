/*
 * 	traqball 2.0
 *	written by Dirk Weber
 *	http://www.eleqtriq.com/
 *	See demo at: http://www.eleqtriq.com/wp-content/static/demos/2011/traqball2011

 *	Copyright (c) 2011 Dirk Weber (http://www.eleqtriq.com)
 *	Licensed under the MIT (http://www.eleqtriq.com/wp-content/uploads/2010/11/mit-license.txt)
 */

(function(){
    var userAgent	= navigator.userAgent.toLowerCase(),
        canTouch 	= "ontouchstart" in window,
        prefix 		= cssPref = "",
        $window   = $(window),
        requestAnimFrame, cancelAnimFrame;

    if(/webkit/gi.test(userAgent)){
        prefix = "-webkit-";
        cssPref = "Webkit";
    }else if(/msie/gi.test(userAgent)){
        prefix = "-ms-";
        cssPref = "ms";
    }else if(/mozilla/gi.test(userAgent)){
        prefix = "-moz-";
        cssPref = "Moz";
    }else if(/opera/gi.test(userAgent)){
        prefix = "-o-";
        cssPref = "O";
    }else{
        prefix = "";
    }

    function bindEvent(target, type, callback, remove){
        //translate events
        var evType 		= type || "touchend",
            mouseEvs 	= ["mousedown", "mouseup", "mousemove"],
            touchEvs 	= ["touchstart", "touchend", "touchmove"],
            remove		= remove || "add";

        evType = canTouch ? evType : mouseEvs[touchEvs.indexOf(type)];

        target[remove+"EventListener"](evType, callback, false);
    }

    function getCoords(eventObj){
        var xTouch,
            yTouch;

        if(eventObj.type.indexOf("mouse") > -1){
            xTouch = eventObj.pageX;
            yTouch = eventObj.pageY;
        }else if(eventObj.type.indexOf("touch") > -1){
            //only do stuff if 1 single finger is used:
            if(eventObj.touches.length === 1){
                var touch	= eventObj.touches[0];
                xTouch		= touch.pageX;
                yTouch		= touch.pageY;
            }
        }

        return [xTouch, yTouch];
    }

    function getStyle(target, prop){
        var style = document.defaultView.getComputedStyle(target, "");
        return style.getPropertyValue(prop);
    }

    requestAnimFrame = (function(){
        return  window[cssPref+"RequestAnimationFrame"] ||
                      function(callback){
                        window.setTimeout(callback, 17);
                      };
        })();

    cancelAnimFrame = (function(){
              return  window[cssPref+"CancelRequestAnimationFrame"] ||
                      clearTimeout;
         })();

    var Traqball = function(confObj){
        this.config = {};
        this.box = null;

        this.setup(confObj);
    };

    Traqball.prototype.disable = function(){
        if(this.box !== null){
            bindEvent(this.box, 'touchstart', this.evHandlers[0], "remove");
            bindEvent(document, 'touchmove', this.evHandlers[1], "remove");
            bindEvent(document, 'touchend', this.evHandlers[2], "remove");
        }
    }

    Traqball.prototype.activate = function(){
        if(this.box !== null){
            bindEvent(this.box, 'touchstart', this.evHandlers[0]);
            bindEvent(document, 'touchmove', this.evHandlers[1], "remove");
            bindEvent(document, 'touchend', this.evHandlers[2], "remove");
        }
    }

    Traqball.prototype.setup = function(conf){
        var THIS			= this,
            radius,					// prepare a variable for storing the radius of our virtual trackball
            stage, 					// the DOM-container of our "rotatable" element
            axis 			= [],	// The rotation-axis
            mouseDownVect 	= [],	// Vector on mousedown
            mouseMoveVect 	= [],	// Vector during mousemove
            startMatrix 	= [],	// Transformation-matrix at the moment of *starting* dragging
            delta 			= 0,
            impulse, pos, w, h, decr, angle, oldAngle, oldTime, curTime;

        (function init (){
            THIS.disable();

            for(var prop in conf){
                THIS.config[prop] = conf[prop];
                }

            stage	= document.getElementById(THIS.config.stage) || document.getElementsByTagname("body")[0];
            pos 	= findPos(stage);
            angle 	= THIS.config.angle || 0;
            impulse	= THIS.config.impulse === false ? false : true;

            // Let's calculate some basic values from "stage" that are necessary for our virtual trackball
            // 1st: determine the radius of our virtual trackball:
            h	= stage.offsetHeight/2,
            w	= stage.offsetWidth/2;

            //take the shortest of both values as radius
            radius = h<w ? h : w;

            //We parse viewport. The first block-element we find will be our "victim" and made rotatable
            for(var i=0, l=stage.childNodes.length; i<l; i++){
                var child = stage.childNodes[i];

                if(child.nodeType === 1){
                    THIS.box = child;
                    break;
                }
            }

            var perspective	= getStyle(stage, prefix+"perspective"),
                pOrigin		= getStyle(stage, prefix+"perspective-origin"),
                bTransform	= getStyle(THIS.box, prefix+"transform");

            //Let's define the start values. If "conf" contains angle or perspective or vector, use them.
            //If not, look for css3d transforms within the CSS.
            //If this fails, let's use some default values.

            if(THIS.config.axis || THIS.config.angle){
                // Normalize the initAxis (initAxis = axis of rotation) because "box" will look distorted if normal is too long
                axis = normalize(THIS.config.axis) || [1,0,0];
                angle = THIS.config.angle || 0;
                // Last but not least we calculate a matrix from the axis and the angle.
                // This matrix will store the initial orientation in 3d-space
                startMatrix = calcMatrix(axis, angle);
            }else if(bTransform !== "none"){
                //already css3d transforms on element?
                startMatrix = bTransform.split(",");

                //Under certain circumstances some browsers report 2d Transforms.
                //Translate them to 3d:
                if(/matrix3d/gi.test(startMatrix[0])){
                    startMatrix[0] = startMatrix[0].replace(/(matrix3d\()/g, "");
                    startMatrix[15] = startMatrix[15].replace(/\)/g, "");
                }else{
                    startMatrix[0] = startMatrix[0].replace(/(matrix\()/g, "");
                    startMatrix[5] = startMatrix[5].replace(/\)/g, "");
                    startMatrix.splice(2,0,0,0);
                    startMatrix.splice(6,0,0,0);
                    startMatrix.splice(8,0,0,0,1,0);
                    startMatrix.splice(14,0,0,1);
                }
                for(var i = 0, l = startMatrix.length; i<l; i++){
                    startMatrix[i] = parseFloat(startMatrix[i]);
                }
            }else{
                axis        = [0,1,0];
                angle       = 0;
                startMatrix = calcMatrix(axis, angle);
            }

            if(THIS.config.perspective){
                stage.style[cssPref+"Perspective"] = THIS.config.perspective;
            }else if(perspective === "none"){
                stage.style[cssPref+"Perspective"] = "700px";
            }

            if(THIS.config.perspectiveOrigin){
                stage.style[cssPref+"PerspectiveOrigin"] = THIS.config.perspectiveOrigin;
            }

            THIS.box.style[cssPref+"Transform"] = "matrix3d("+ startMatrix+")";
            bindEvent(THIS.box, 'touchstart', startrotation);

            THIS.evHandlers = [startrotation, rotate, finishrotation];
        })();


        function startrotation(e){
            if(delta !== 0){stopSlide();};
            e.preventDefault();

            mouseDownVect 	= calcZvector(getCoords(e));
            oldTime			= curTime = new Date().getTime();
            oldAngle 		= angle;

            bindEvent(THIS.box,'touchstart', startrotation, "remove");
            bindEvent(document, 'touchmove', rotate);
            bindEvent(document, 'touchend', finishrotation);
        }

        function finishrotation(e){
            var stopMatrix;

            bindEvent(document, 'touchmove', rotate, "remove");
            bindEvent(document, 'touchend', finishrotation, "remove");
            bindEvent(THIS.box, 'touchstart', startrotation);
            calcSpeed();
            if( impulse && delta > 0){
                requestAnimFrame(slide);
            }else if(!(isNaN(axis[0]) || isNaN(axis[1]) || isNaN(axis[2]))) {
                stopSlide();
            }
        }

        function cleanupMatrix(){
            // Clean up when finishing rotation. Only thing to do: create a new "initial" matrix for the next rotation.
            // If we don't, the object will flip back to the position at launch every time the user starts dragging.
            // Therefore we must:
            // 1. calculate a matrix from axis and the current angle
            // 2. Create a new startmatrix by combining current startmatrix and stopmatrix to a new matrix.
            // Matrices can be combined by multiplication, so what are we waiting for?
            stopMatrix	= calcMatrix(axis, angle);
            startMatrix	= multiplyMatrix(startMatrix,stopMatrix);
        }

        // The rotation:
        function rotate(e){
            var eCoords	= getCoords(e);
            e.preventDefault();

            oldTime = curTime;
            oldAngle = angle;

            // Calculate the currrent z-component of the 3d-vector on the virtual trackball
            mouseMoveVect = calcZvector(eCoords);

            // We already calculated the z-vector-component on mousedown and the z-vector-component during mouse-movement.
            // We will use them to retrieve the current rotation-axis
            // (the normal-vector perpendiular to mouseDownVect and mouseMoveVect).
            axis[0] = mouseDownVect[1]*mouseMoveVect[2]-mouseDownVect[2]*mouseMoveVect[1];
            axis[1] = mouseDownVect[2]*mouseMoveVect[0]-mouseDownVect[0]*mouseMoveVect[2];
            axis[2] = mouseDownVect[0]*mouseMoveVect[1]-mouseDownVect[1]*mouseMoveVect[0];
            axis	= normalize(axis);

            // Now that we have the normal, we need the angle of the rotation.
            // Easy to find by calculating the angle between mouseDownVect and mouseMoveVect:
            angle = calcAngle(mouseDownVect, mouseMoveVect);

            //Only one thing left to do: Update the position of the box by applying a new transform:
            // 2 transforms will be applied: the current rotation 3d and the start-matrix
            THIS.box.style[cssPref+"Transform"] = "rotate3d("+ axis+","+angle+"rad) matrix3d("+startMatrix+")";

            $window.trigger('traqBallRotate');

            curTime = new Date().getTime();
        }

        function calcSpeed(){
            var dw 	= angle - oldAngle;
                dt 	= curTime-oldTime;

            delta 	= Math.abs(dw*17/dt);

            if(isNaN(delta)){
                delta = 0;
            }else if(delta > 0.2){
                delta = 0.2;
            }
        }

        function slide(){
            angle+= delta;
            decr = 0.01*Math.sqrt(delta);
            delta = delta > 0 ? delta-decr : 0;

            THIS.box.style[cssPref+"Transform"] = "rotate3d("+ axis+","+angle+"rad) matrix3d("+startMatrix+")";
            $window.trigger('traqBallRotate');

            if (delta === 0){
                stopSlide();
            }else{
                requestAnimFrame(slide);
            }
        }

        function stopSlide(){
            cancelAnimFrame(slide);
            cleanupMatrix();
            oldAngle = angle = 0;
            delta = 0;
        }

        //Some stupid matrix-multiplication.
        function multiplyMatrix(m1, m2){
            var matrix = [];

            matrix[0]	= m1[0]*m2[0]+m1[1]*m2[4]+m1[2]*m2[8]+m1[3]*m2[12];
            matrix[1]	= m1[0]*m2[1]+m1[1]*m2[5]+m1[2]*m2[9]+m1[3]*m2[13];
            matrix[2]	= m1[0]*m2[2]+m1[1]*m2[6]+m1[2]*m2[10]+m1[3]*m2[14];
            matrix[3]	= m1[0]*m2[3]+m1[1]*m2[7]+m1[2]*m2[11]+m1[3]*m2[15];
            matrix[4]	= m1[4]*m2[0]+m1[5]*m2[4]+m1[6]*m2[8]+m1[7]*m2[12];
            matrix[5]	= m1[4]*m2[1]+m1[5]*m2[5]+m1[6]*m2[9]+m1[7]*m2[13];
            matrix[6]	= m1[4]*m2[2]+m1[5]*m2[6]+m1[6]*m2[10]+m1[7]*m2[14];
            matrix[7]	= m1[4]*m2[3]+m1[5]*m2[7]+m1[6]*m2[11]+m1[7]*m2[15];
            matrix[8]	= m1[8]*m2[0]+m1[9]*m2[4]+m1[10]*m2[8]+m1[11]*m2[12];
            matrix[9]	= m1[8]*m2[1]+m1[9]*m2[5]+m1[10]*m2[9]+m1[11]*m2[13];
            matrix[10]	= m1[8]*m2[2]+m1[9]*m2[6]+m1[10]*m2[10]+m1[11]*m2[14];
            matrix[11]	= m1[8]*m2[3]+m1[9]*m2[7]+m1[10]*m2[11]+m1[11]*m2[15];
            matrix[12]	= m1[12]*m2[0]+m1[13]*m2[4]+m1[14]*m2[8]+m1[15]*m2[12];
            matrix[13]	= m1[12]*m2[1]+m1[13]*m2[5]+m1[14]*m2[9]+m1[15]*m2[13];
            matrix[14]	= m1[12]*m2[2]+m1[13]*m2[6]+m1[14]*m2[10]+m1[15]*m2[14];
            matrix[15]	= m1[12]*m2[3]+m1[13]*m2[7]+m1[14]*m2[11]+m1[15]*m2[15];

            return matrix;
        }

        // This function will calculate a z-component for our 3D-vector from the mouse x and y-coordinates
        // (the corresponding point on our virtual trackball):
        function calcZvector(coords){
            var x 		= coords[0] - pos[0],
                y 		= coords[1] - pos[1],
                vector 	= [(x/radius-1), (y/radius-1)],
                z 		= 1 - vector[0]*vector[0] - vector[1]*vector[1];

             // Make sure that dragging stops when z gets a negative value:
            vector[2] 	= z > 0 ? Math.sqrt(z) : 0;

            return vector;
        }

        // Normalization recalculates all coordinates in a way that the resulting vector has a length of "1".
        // We achieve this by dividing the x, y and z-coordinates by the vector's length
        function normalize(vect){
            var length = Math.sqrt( vect[0]*vect[0] + vect[1]*vect[1] + vect[2]*vect[2] );

            vect[0]/= length;
            vect[1]/= length;
            vect[2]/= length;

            return vect;
        }

        // Calculate the angle between 2 vectors.
        function calcAngle(vect_1, vect_2){
            var numerator 	= 	vect_1[0]*vect_2[0] + vect_1[1]*vect_2[1] + vect_1[2]*vect_2[2],
                denominator	= 	Math.sqrt(vect_1[0]*vect_1[0] + vect_1[1]*vect_1[1] + vect_1[2]*vect_1[2])*
                                Math.sqrt(vect_2[0]*vect_2[0] + vect_2[1]*vect_2[1] + vect_2[2]*vect_2[2]),
                angle		=	Math.acos(numerator/denominator);

            return angle;
        }

        function calcMatrix(vector, angle){
            // calculate transformation-matrix from a vector[x,y,z] and an angle
            var x		= vector[0],
                y		= vector[1],
                z		= vector[2],
                sin		= Math.sin(angle),
                cos		= Math.cos(angle),
                cosmin	= 1-cos,
                matrix	= [(cos+x*x*cosmin), (y*x*cosmin+z*sin),(z*x*cosmin-y*sin),0,
                          (x*y*cosmin-z*sin), (cos+y*y*cosmin), (z*y*cosmin+x*sin),0,
                          (x*z*cosmin+y*sin), (y*z*cosmin-x*sin), (cos+z*z*cosmin),0,
                          0,0,0,1];

            return matrix;
        }

        //findPos-script by www.quirksmode.org
        function findPos(obj) {
            var curleft = 0,
                curtop 	= 0;

            if (obj.offsetParent) {
                do {
                    curleft += obj.offsetLeft;
                    curtop += obj.offsetTop;
                } while (obj = obj.offsetParent);

                return [curleft,curtop];
            }
        }
    }

    window.Traqball = Traqball;
})();


/*
 * Photon
 * http://photon.attasi.com
 *
 * Licensed under the MIT license.
 * Copyright 2012 Tom Giannattasio
 */

var Photon={version:"0.0.3",degToRad:function(a){return a*Math.PI/180},radToDeg:function(a){return a*180/Math.PI},getRotationVector:function(b,a){var e=b.rotate(a.x,Line.create([0,0,0],[1,0,0]));var c=e.rotate(a.y,Line.create([0,0,0],[0,1,0]));var d=c.rotate(a.z,Line.create([0,0,0],[0,0,1]));return d},getTransformString:function(){if(Photon.transformString){return Photon.transformString}var c;var d=["transform","webkitTransform","MozTransform","msTransform","OTransform"];var b=document.createElement("div");for(var a=0;a<d.length;a++){if(b.style[d[a]]==""){c=d[a]}}Photon.transformString=c;return c},buildMatrix:function(b){var a=new FirminCSSMatrix(b);a.m11=a.m11*10000000000000000;a.m12=a.m12*10000000000000000;a.m13=a.m13*10000000000000000;a.m14=a.m14*10000000000000000;a.m21=a.m21*10000000000000000;a.m22=a.m22*10000000000000000;a.m23=a.m23*10000000000000000;a.m24=a.m24*10000000000000000;a.m31=a.m31*10000000000000000;a.m32=a.m32*10000000000000000;a.m33=a.m33*10000000000000000;a.m34=a.m34*10000000000000000;a.m41=a.m41*10000000000000000;a.m42=a.m42*10000000000000000;a.m43=a.m43*10000000000000000;a.m44=a.m44*10000000000000000;return a}};Photon.Light=function(c,b,a){this.moveTo(c||0,b||0,a||100);this.calculateVector()};Photon.Light.prototype={moveTo:function(a,c,b){this.x=a;this.y=c;this.z=b;this.calculateVector()},calculateVector:function(){this.magnitude=Math.sqrt((this.x*this.x)+(this.y*this.y)+(this.z*this.z));this.vector=$V([this.x/this.magnitude,this.y/this.magnitude,this.z/this.magnitude])}};Photon.Face=function(d,b,a,c){this.element=d;this.maxShade=b||0.5;this.maxTint=a||0;this.isBackfaced=c||false;this.shaderElement=new Photon.ShaderElement(this.element);this.element.insertBefore(this.shaderElement,this.element.firstChild);this.transformString=Photon.getTransformString();this.getRotations()};Photon.Face.prototype={getRotations:function(){var b=window.getComputedStyle(this.element)[this.transformString]||"matrix3d(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)";this.matrix=Photon.buildMatrix(b);var a=this.matrix.decompose();this.rotations={x:a.rotate.x,y:a.rotate.y,z:a.rotate.z};this.vector=Photon.getRotationVector($V([0,0,1]),this.rotations)},render:function(b,h,g){if(h){this.getRotations()}var e;if(g){e=Photon.getRotationVector(this.vector,g)}else{e=this.vector}this.angleFrom=Photon.radToDeg(b.vector.angleFrom(e));var f;var d=this.isBackfaced?this.angleFrom/180:this.angleFrom/90;if(this.isBackfaced&&d>0.5){d=1-d}var c=Math.abs(this.maxShade+this.maxTint);var a=c*d;this.rangedPercentage=a;if(a<=this.maxTint){f="rgba(255, 255, 255, "+Math.abs(this.maxTint-a)+")"}else{f="rgba(0, 0, 0, "+Math.abs(a-this.maxTint)+")"}this.shaderElement.style.background=f},setMaxShade:function(a){this.maxShade=a},setMaxTint:function(a){this.maxTint=a}};Photon.ShaderElement=function(a){var b=document.createElement("div");b.className="photon-shader";b.style.position="absolute";b.style.top="0";b.style.left="0";b.style.width=window.getComputedStyle(a).width;b.style.height=window.getComputedStyle(a).height;return b};Photon.FaceGroup=function(f,a,c,b,d){this.element=f;this.faces=[];this.transformString=Photon.getTransformString();var g=a;for(var e=0;e<g.length;e++){this.faces[e]=new Photon.Face(g[e],c,b,d)}};Photon.FaceGroup.prototype={getRotations:function(){var b=window.getComputedStyle(this.element)[this.transformString]||"matrix3d(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)";this.matrix=Photon.buildMatrix(b);var a=this.matrix.decompose();this.rotations={x:a.rotate.x,y:a.rotate.y,z:a.rotate.z};this.vector=Photon.getRotationVector($V([0,0,1]),this.rotations)},render:function(a,d,c){if(d){this.getRotations()}this.angleFrom=Photon.radToDeg(a.vector.angleFrom(this.vector));for(var b=0,e=this.faces.length;b<e;b++){this.faces[b].render(a,c,this.rotations)}},setMaxShade:function(b){for(var a=0;a<this.faces.length;a++){this.faces[a].setMaxShade(b)}},setMaxTint:function(b){for(var a=0;a<this.faces.length;a++){this.faces[a].setMaxTint(b)}}};var Sylvester={version:"0.1.3",precision:0.000001};function Vector(){}Vector.prototype={modulus:function(){return Math.sqrt(this.dot(this))},dup:function(){return Vector.create(this.elements)},each:function(c){var d=this.elements.length,a=d,b;do{b=a-d;c(this.elements[b],b+1)}while(--d)},angleFrom:function(d){var e=d.elements||d;var c=this.elements.length,f=c,g;if(c!=e.length){return null}var a=0,l=0,h=0;this.each(function(k,m){a+=k*e[m-1];l+=k*k;h+=e[m-1]*e[m-1]});l=Math.sqrt(l);h=Math.sqrt(h);if(l*h===0){return null}var b=a/(l*h);if(b<-1){b=-1}if(b>1){b=1}return Math.acos(b)},dot:function(b){var a=b.elements||b;var c,d=0,e=this.elements.length;if(e!=a.length){return null}do{d+=this.elements[e-1]*a[e-1]}while(--e);return d},rotate:function(c,e){var b,d,a,h,g;switch(this.elements.length){case 2:b=e.elements||e;if(b.length!=2){return null}d=Matrix.Rotation(c).elements;a=this.elements[0]-b[0];h=this.elements[1]-b[1];return Vector.create([b[0]+d[0][0]*a+d[0][1]*h,b[1]+d[1][0]*a+d[1][1]*h]);break;case 3:if(!e.direction){return null}var f=e.pointClosestTo(this).elements;d=Matrix.Rotation(c,e.direction).elements;a=this.elements[0]-f[0];h=this.elements[1]-f[1];g=this.elements[2]-f[2];return Vector.create([f[0]+d[0][0]*a+d[0][1]*h+d[0][2]*g,f[1]+d[1][0]*a+d[1][1]*h+d[1][2]*g,f[2]+d[2][0]*a+d[2][1]*h+d[2][2]*g]);break;default:return null}},setElements:function(a){this.elements=(a.elements||a).slice();return this}};Vector.create=function(b){var a=new Vector();return a.setElements(b)};var $V=Vector.create;function Line(){}Line.prototype={distanceFrom:function(e){if(e.normal){return e.distanceFrom(this)}if(e.direction){if(this.isParallelTo(e)){return this.distanceFrom(e.anchor)}var k=this.direction.cross(e.direction).toUnitVector().elements;var c=this.anchor.elements,b=e.anchor.elements;return Math.abs((c[0]-b[0])*k[0]+(c[1]-b[1])*k[1]+(c[2]-b[2])*k[2])}else{var f=e.elements||e;var c=this.anchor.elements,a=this.direction.elements;var n=f[0]-c[0],l=f[1]-c[1],g=(f[2]||0)-c[2];var m=Math.sqrt(n*n+l*l+g*g);if(m===0){return 0}var h=(n*a[0]+l*a[1]+g*a[2])/m;var d=1-h*h;return Math.abs(m*Math.sqrt(d<0?0:d))}},contains:function(a){var b=this.distanceFrom(a);return(b!==null&&b<=Sylvester.precision)},pointClosestTo:function(s){if(s.direction){if(this.intersects(s)){return this.intersectionWith(s)}if(this.isParallelTo(s)){return null}var u=this.direction.elements,t=s.direction.elements;var f=u[0],e=u[1],c=u[2],q=t[0],o=t[1],m=t[2];var r=(c*q-f*m),p=(f*o-e*q),n=(e*m-c*o);var l=Vector.create([r*m-p*o,p*q-n*m,n*o-r*q]);var h=Plane.create(s.anchor,l);return h.intersectionWith(this)}else{var h=s.elements||s;if(this.contains(h)){return Vector.create(h)}var v=this.anchor.elements,u=this.direction.elements;var f=u[0],e=u[1],c=u[2],d=v[0],b=v[1],a=v[2];var r=f*(h[1]-b)-e*(h[0]-d),p=e*((h[2]||0)-a)-c*(h[1]-b),n=c*(h[0]-d)-f*((h[2]||0)-a);var g=Vector.create([e*r-c*n,c*p-f*r,f*n-e*p]);var w=this.distanceFrom(h)/g.modulus();return Vector.create([h[0]+g.elements[0]*w,h[1]+g.elements[1]*w,(h[2]||0)+g.elements[2]*w])}},rotate:function(p,q){if(typeof(q.direction)=="undefined"){q=Line.create(q.to3D(),Vector.k)}var g=Matrix.Rotation(p,q.direction).elements;var b=q.pointClosestTo(this.anchor).elements;var d=this.anchor.elements,a=this.direction.elements;var l=b[0],k=b[1],h=b[2],f=d[0],e=d[1],c=d[2];var o=f-l,n=e-k,m=c-h;return Line.create([l+g[0][0]*o+g[0][1]*n+g[0][2]*m,k+g[1][0]*o+g[1][1]*n+g[1][2]*m,h+g[2][0]*o+g[2][1]*n+g[2][2]*m],[g[0][0]*a[0]+g[0][1]*a[1]+g[0][2]*a[2],g[1][0]*a[0]+g[1][1]*a[1]+g[1][2]*a[2],g[2][0]*a[0]+g[2][1]*a[1]+g[2][2]*a[2]])},setVectors:function(a,c){a=Vector.create(a);c=Vector.create(c);if(a.elements.length==2){a.elements.push(0)}if(c.elements.length==2){c.elements.push(0)}if(a.elements.length>3||c.elements.length>3){return null}var b=c.modulus();if(b===0){return null}this.anchor=a;this.direction=Vector.create([c.elements[0]/b,c.elements[1]/b,c.elements[2]/b]);return this}};Line.create=function(b,c){var a=new Line();return a.setVectors(b,c)};function Matrix(){}Matrix.prototype={setElements:function(h){var m,a=h.elements||h;if(typeof(a[0][0])!="undefined"){var d=a.length,f=d,b,c,l;this.elements=[];do{m=f-d;b=a[m].length;c=b;this.elements[m]=[];do{l=c-b;this.elements[m][l]=a[m][l]}while(--b)}while(--d);return this}var e=a.length,g=e;this.elements=[];do{m=g-e;this.elements.push([a[m]])}while(--e);return this}};Matrix.create=function(a){var b=new Matrix();return b.setElements(a)};Matrix.Rotation=function(b,k){if(!k){return Matrix.create([[Math.cos(b),-Math.sin(b)],[Math.sin(b),Math.cos(b)]])}var d=k.dup();if(d.elements.length!=3){return null}var h=d.modulus();var l=d.elements[0]/h,g=d.elements[1]/h,f=d.elements[2]/h;var n=Math.sin(b),e=Math.cos(b),m=1-e;return Matrix.create([[m*l*l+e,m*l*g-n*f,m*l*f+n*g],[m*l*g+n*f,m*g*g+e,m*g*f-n*l],[m*l*f-n*g,m*g*f+n*l,m*f*f+e]])};FirminCSSMatrix=function(a){this.m11=this.m22=this.m33=this.m44=1;this.m12=this.m13=this.m14=this.m21=this.m23=this.m24=this.m31=this.m32=this.m34=this.m41=this.m42=this.m43=0;if(typeof a=="string"){this.setMatrixValue(a)}};FirminCSSMatrix.displayName="FirminCSSMatrix";FirminCSSMatrix.degreesToRadians=function(a){return a*Math.PI/180};FirminCSSMatrix.prototype.isAffine=function(){return this.m13===0&&this.m14===0&&this.m23===0&&this.m24===0&&this.m31===0&&this.m32===0&&this.m33===1&&this.m34===0&&this.m43===0&&this.m44===1};FirminCSSMatrix.prototype.setMatrixValue=function(g){g=g.trim();var b=g.match(/^matrix(3d)?\(\s*(.+)\s*\)$/),f,h,a,e,d,c;if(!b){return}f=!!b[1];h=b[2].split(/\s*,\s*/);a=h.length;e=new Array(a);if((f&&a!==16)||!(f||a===6)){return}for(d=0;d<a;d++){c=h[d];if(c.match(/^-?\d+(\.\d+)?$/)){e[d]=parseFloat(c)}else{return}}for(d=0;d<a;d++){point=f?("m"+(Math.floor(d/4)+1))+(d%4+1):String.fromCharCode(d+97);this[point]=e[d]}};FirminCSSMatrix.prototype.toString=function(){var a=this,b,c;if(this.isAffine()){c="matrix(";b=["a","b","c","d","e","f"]}else{c="matrix3d(";b=["m11","m12","m13","m14","m21","m22","m23","m24","m31","m32","m33","m34","m41","m42","m43","m44"]}return c+b.map(function(d){return a[d].toFixed(6)}).join(", ")+")"};var CSSMatrixDecomposed=function(c){c===undefined?c={}:null;var b={perspective:null,translate:null,skew:null,scale:null,rotate:null};for(var a in b){this[a]=c[a]?c[a]:new Vector4()}this.tween=function(d,f,k){if(k===undefined){k=function(e){return e}}if(!d){d=new CSSMatrixDecomposed(new FirminCSSMatrix().decompose())}var l=new CSSMatrixDecomposed(),h=index=null,g="";f=k(f);for(index in b){for(h in {x:"x",y:"y",z:"z",w:"w"}){l[index][h]=(this[index][h]+(d[index][h]-this[index][h])*f).toFixed(5)}}g="matrix3d(1,0,0,0, 0,1,0,0, 0,0,1,0, "+l.perspective.x+", "+l.perspective.y+", "+l.perspective.z+", "+l.perspective.w+") translate3d("+l.translate.x+"px, "+l.translate.y+"px, "+l.translate.y+"px) rotateX("+l.rotate.x+"rad) rotateY("+l.rotate.y+"rad) rotateZ("+l.rotate.z+"rad) matrix3d(1,0,0,0, 0,1,0,0, 0,"+l.skew.z+",1,0, 0,0,0,1) matrix3d(1,0,0,0, 0,1,0,0, "+l.skew.y+",0,1,0, 0,0,0,1) matrix3d(1,0,0,0, "+l.skew.x+",1,0,0, 0,0,1,0, 0,0,0,1) scale3d("+l.scale.x+", "+l.scale.y+", "+l.scale.z+")";try{l=new FirminCSSMatrix(g);return l}catch(m){console.error("Invalid matrix string: "+g);return""}}};var Vector4=function(a,d,c,b){this.x=a?a:0;this.y=d?d:0;this.z=c?c:0;this.w=b?b:0;this.checkValues=function(){this.x=this.x?this.x:0;this.y=this.y?this.y:0;this.z=this.z?this.z:0;this.w=this.w?this.w:0};this.length=function(){this.checkValues();return Math.sqrt(this.x*this.x+this.y*this.y+this.z*this.z)};this.normalise=function(){var e=this.length(),f=new Vector4(this.x/e,this.y/e,this.z/e);return f};this.dot=function(e){return this.x*e.x+this.y*e.y+this.z*e.z+this.w*e.w};this.cross=function(e){return new Vector4(this.y*e.z-this.z*e.y,this.z*e.x-this.x*e.z,this.x*e.y-this.y*e.x)};this.combine=function(g,f,e){return new Vector4((f*this.x)+(e*g.x),(f*this.y)+(e*g.y),(f*this.z)+(e*g.z))}};FirminCSSMatrix.prototype.determinant=function(){return this.m14*this.m23*this.m32*this.m41-this.m13*this.m24*this.m32*this.m41-this.m14*this.m22*this.m33*this.m41+this.m12*this.m24*this.m33*this.m41+this.m13*this.m22*this.m34*this.m41-this.m12*this.m23*this.m34*this.m41-this.m14*this.m23*this.m31*this.m42+this.m13*this.m24*this.m31*this.m42+this.m14*this.m21*this.m33*this.m42-this.m11*this.m24*this.m33*this.m42-this.m13*this.m21*this.m34*this.m42+this.m11*this.m23*this.m34*this.m42+this.m14*this.m22*this.m31*this.m43-this.m12*this.m24*this.m31*this.m43-this.m14*this.m21*this.m32*this.m43+this.m11*this.m24*this.m32*this.m43+this.m12*this.m21*this.m34*this.m43-this.m11*this.m22*this.m34*this.m43-this.m13*this.m22*this.m31*this.m44+this.m12*this.m23*this.m31*this.m44+this.m13*this.m21*this.m32*this.m44-this.m11*this.m23*this.m32*this.m44-this.m12*this.m21*this.m33*this.m44+this.m11*this.m22*this.m33*this.m44};FirminCSSMatrix.prototype.decompose=function(){var a=new FirminCSSMatrix(this.toString()),b=rightHandSide=inversePerspectiveMatrix=transposedInversePerspectiveMatrix=perspective=translate=row=i=scale=skew=pdum3=rotate=null;if(a.m33==0){return new CSSMatrixDecomposed(new FirminCSSMatrix().decompose())}for(i=1;i<=4;i++){for(j=1;j<=4;j++){a["m"+i+j]/=a.m44}}b=a;for(i=1;i<=3;i++){b["m"+i+"4"]=0}b.m44=1;if(b.determinant()==0){return new CSSMatrixDecomposed(new FirminCSSMatrix().decompose())}if(a.m14!=0||a.m24!=0||a.m34!=0){rightHandSide=new Vector4(a.m14,a.m24,a.m34,a.m44);inversePerspectiveMatrix=b.inverse();transposedInversePerspectiveMatrix=inversePerspectiveMatrix.transpose();perspective=transposedInversePerspectiveMatrix.transformVector(rightHandSide);a.m14=a.m24=a.m34=0;a.m44=1}else{perspective=new Vector4(0,0,0,1)}translate=new Vector4(a.m41,a.m42,a.m43);a.m41=0;a.m42=0;a.m43=0;row=[new Vector4(),new Vector4(),new Vector4()];for(i=1;i<=3;i++){row[i-1].x=a["m"+i+"1"];row[i-1].y=a["m"+i+"2"];row[i-1].z=a["m"+i+"3"]}scale=new Vector4();skew=new Vector4();scale.x=row[0].length();row[0]=row[0].normalise();skew.x=row[0].dot(row[1]);row[1]=row[1].combine(row[0],1,-skew.x);scale.y=row[1].length();row[1]=row[1].normalise();skew.x/=scale.y;skew.y=row[0].dot(row[2]);row[2]=row[2].combine(row[0],1,-skew.y);skew.z=row[1].dot(row[2]);row[2]=row[2].combine(row[1],1,-skew.z);scale.z=row[2].length();row[2]=row[2].normalise();skew.y/=scale.z;skew.y/=scale.z;pdum3=row[1].cross(row[2]);if(row[0].dot(pdum3)<0){for(i=0;i<3;i++){scale.x*=-1;row[i].x*=-1;row[i].y*=-1;row[i].z*=-1}}rotate=new Vector4();rotate.y=Math.asin(-row[0].z);if(Math.cos(rotate.y)!=0){rotate.x=Math.atan2(row[1].z,row[2].z);rotate.z=Math.atan2(row[0].y,row[0].x)}else{rotate.x=Math.atan2(-row[2].x,row[1].y);rotate.z=0}return new CSSMatrixDecomposed({perspective:perspective,translate:translate,skew:skew,scale:scale,rotate:rotate})};
