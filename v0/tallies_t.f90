module xml_data_tallies_t
   use READ_XML_PRIMITIVES
   use WRITE_XML_PRIMITIVES
   use XMLPARSE
   implicit none
   save
   integer, private :: lurep_
   logical, private :: strict_

type mesh_xml
   integer                                         :: id
   character(len=12)                                :: type
   integer, dimension(:), pointer                  :: dimension => null()
   real(kind=kind(1.0d0)), dimension(:), pointer   :: lower_left => null()
   real(kind=kind(1.0d0)), dimension(:), pointer   :: upper_right => null()
   real(kind=kind(1.0d0)), dimension(:), pointer   :: width => null()
end type mesh_xml

type filter_xml
   character(len=20)                                :: type
   character(len=20), dimension(:), pointer         :: bins => null()
   character(len=20)                                :: groups
end type filter_xml

type tally_xml
   integer                                         :: id
   character(len=52)                                :: label
   character(len=12)                                :: estimator
   type(filter_xml), dimension(:), pointer         :: filter => null()
   character(len=20), dimension(:), pointer         :: scores => null()
   character(len=12), dimension(:), pointer         :: nuclides => null()
   type(filter_xml), dimension(:), pointer         :: filters => null()
end type tally_xml
   type(mesh_xml), dimension(:), pointer           :: mesh_ => null()
   type(tally_xml), dimension(:), pointer          :: tally_ => null()
   character(len=3)                                :: separate_
contains
subroutine read_xml_type_mesh_xml_array( &
      info, tag, endtag, attribs, noattribs, data, nodata, &
      dvar, has_dvar )
   type(XML_PARSE)                                 :: info
   character(len=*), intent(inout)                 :: tag
   logical, intent(inout)                          :: endtag
   character(len=*), dimension(:,:), intent(inout) :: attribs
   integer, intent(inout)                          :: noattribs
   character(len=*), dimension(:), intent(inout)   :: data
   integer, intent(inout)                          :: nodata
   type(mesh_xml), dimension(:), pointer :: dvar
   logical, intent(inout)                       :: has_dvar

   integer                                      :: newsize
   type(mesh_xml), dimension(:), pointer :: newvar

   newsize = size(dvar) + 1
   allocate( newvar(1:newsize) )
   newvar(1:newsize-1) = dvar
   deallocate( dvar )
   dvar => newvar

   call read_xml_type_mesh_xml( info, tag, endtag, attribs, noattribs, data, nodata, &
              dvar(newsize), has_dvar )
end subroutine read_xml_type_mesh_xml_array

subroutine read_xml_type_mesh_xml( info, starttag, endtag, attribs, noattribs, data, nodata, &
              dvar, has_dvar )
   type(XML_PARSE)                                 :: info
   character(len=*), intent(in)                    :: starttag
   logical, intent(inout)                          :: endtag
   character(len=*), dimension(:,:), intent(inout) :: attribs
   integer, intent(inout)                          :: noattribs
   character(len=*), dimension(:), intent(inout)   :: data
   integer, intent(inout)                          :: nodata
   type(mesh_xml), intent(inout)  :: dvar
   logical, intent(inout)                       :: has_dvar

   integer                                      :: att_
   integer                                      :: noatt_
   logical                                      :: error
   logical                                      :: endtag_org
   character(len=len(starttag))                 :: tag
   logical                                         :: has_id
   logical                                         :: has_type
   logical                                         :: has_dimension
   logical                                         :: has_lower_left
   logical                                         :: has_upper_right
   logical                                         :: has_width
   has_id                               = .false.
   has_type                             = .false.
   has_dimension                        = .false.
   has_lower_left                       = .false.
   has_upper_right                      = .false.
   has_width                            = .false.
   call init_xml_type_mesh_xml(dvar)
   has_dvar = .true.
   error  = .false.
   att_   = 0
   noatt_ = noattribs+1
   endtag_org = endtag
   do
      if ( nodata /= 0 ) then
         noattribs = 0
         tag = starttag
      elseif ( att_ < noatt_ .and. noatt_ > 1 ) then
         att_      = att_ + 1
         if ( att_ <= noatt_-1 ) then
            tag       = attribs(1,att_)
            data(1)   = attribs(2,att_)
            noattribs = 0
            nodata    = 1
            endtag    = .false.
         else
            tag       = starttag
            noattribs = 0
            nodata    = 0
            endtag    = .true.
            cycle
         endif
      else
         if ( endtag_org ) then
            return
         else
            call xml_get( info, tag, endtag, attribs, noattribs, data, nodata )
            if ( xml_error(info) ) then
               write(lurep_,*) 'Error reading input file!'
               error = .true.
               return
            endif
         endif
      endif
      if ( endtag .and. tag == starttag ) then
         exit
      endif
      if ( endtag .and. noattribs == 0 ) then
         if ( xml_ok(info) ) then
            cycle
         else
            exit
         endif
      endif
      select case( tag )
      case('id')
         call read_xml_integer( &
            info, tag, endtag, attribs, noattribs, data, nodata, &
            dvar%id, has_id )
      case('type')
         call read_xml_word( &
            info, tag, endtag, attribs, noattribs, data, nodata, &
            dvar%type, has_type )
      case('dimension')
         call read_xml_integer_array( &
            info, tag, endtag, attribs, noattribs, data, nodata, &
            dvar%dimension, has_dimension )
      case('lower_left')
         call read_xml_double_array( &
            info, tag, endtag, attribs, noattribs, data, nodata, &
            dvar%lower_left, has_lower_left )
      case('upper_right')
         call read_xml_double_array( &
            info, tag, endtag, attribs, noattribs, data, nodata, &
            dvar%upper_right, has_upper_right )
      case('width')
         call read_xml_double_array( &
            info, tag, endtag, attribs, noattribs, data, nodata, &
            dvar%width, has_width )
      case ('comment', '!--')
         ! Simply ignore
      case default
         if ( strict_ ) then
            error = .true.
            call xml_report_errors( info, &
               'Unknown or wrongly placed tag: ' // trim(tag))
         endif
      end select
      nodata = 0
      if ( .not. xml_ok(info) ) exit
   end do
   if ( .not. has_id ) then
      has_dvar = .false.
      call xml_report_errors(info, 'Missing data on id')
   endif
   if ( .not. has_type ) then
      has_dvar = .false.
      call xml_report_errors(info, 'Missing data on type')
   endif
   if ( .not. has_dimension ) then
      has_dvar = .false.
      call xml_report_errors(info, 'Missing data on dimension')
   endif
   if ( .not. has_lower_left ) then
      has_dvar = .false.
      call xml_report_errors(info, 'Missing data on lower_left')
   endif
   if ( .not. has_upper_right ) then
      has_dvar = .false.
      call xml_report_errors(info, 'Missing data on upper_right')
   endif
   if ( .not. has_width ) then
      has_dvar = .false.
      call xml_report_errors(info, 'Missing data on width')
   endif
end subroutine read_xml_type_mesh_xml
subroutine init_xml_type_mesh_xml_array( dvar )
   type(mesh_xml), dimension(:), pointer :: dvar
   if ( associated( dvar ) ) then
      deallocate( dvar )
   endif
   allocate( dvar(0) )
end subroutine init_xml_type_mesh_xml_array
subroutine init_xml_type_mesh_xml(dvar)
   type(mesh_xml) :: dvar
end subroutine init_xml_type_mesh_xml
subroutine write_xml_type_mesh_xml_array( &
      info, tag, indent, dvar )
   type(XML_PARSE)                                 :: info
   character(len=*), intent(in)                    :: tag
   integer                                         :: indent
   type(mesh_xml), dimension(:)        :: dvar
   integer                                         :: i
   do i = 1,size(dvar)
       call write_xml_type_mesh_xml( info, tag, indent, dvar(i) )
   enddo
end subroutine write_xml_type_mesh_xml_array

subroutine write_xml_type_mesh_xml( &
      info, tag, indent, dvar )
   type(XML_PARSE)                                 :: info
   character(len=*), intent(in)                    :: tag
   integer                                         :: indent
   type(mesh_xml)                      :: dvar
   character(len=100)                              :: indentation
   indentation = ' '
   write(info%lun, '(4a)' ) indentation(1:min(indent,100)),&
       '<',trim(tag), '>'
   call write_to_xml_integer( info, 'id', indent+3, dvar%id)
   call write_to_xml_word( info, 'type', indent+3, dvar%type)
   call write_to_xml_integer_array( info, 'dimension', indent+3, dvar%dimension)
   call write_to_xml_double_array( info, 'lower_left', indent+3, dvar%lower_left)
   call write_to_xml_double_array( info, 'upper_right', indent+3, dvar%upper_right)
   call write_to_xml_double_array( info, 'width', indent+3, dvar%width)
   write(info%lun,'(4a)') indentation(1:min(indent,100)), &
       '</' //trim(tag) // '>'
end subroutine write_xml_type_mesh_xml

subroutine read_xml_type_filter_xml_array( &
      info, tag, endtag, attribs, noattribs, data, nodata, &
      dvar, has_dvar )
   type(XML_PARSE)                                 :: info
   character(len=*), intent(inout)                 :: tag
   logical, intent(inout)                          :: endtag
   character(len=*), dimension(:,:), intent(inout) :: attribs
   integer, intent(inout)                          :: noattribs
   character(len=*), dimension(:), intent(inout)   :: data
   integer, intent(inout)                          :: nodata
   type(filter_xml), dimension(:), pointer :: dvar
   logical, intent(inout)                       :: has_dvar

   integer                                      :: newsize
   type(filter_xml), dimension(:), pointer :: newvar

   newsize = size(dvar) + 1
   allocate( newvar(1:newsize) )
   newvar(1:newsize-1) = dvar
   deallocate( dvar )
   dvar => newvar

   call read_xml_type_filter_xml( info, tag, endtag, attribs, noattribs, data, nodata, &
              dvar(newsize), has_dvar )
end subroutine read_xml_type_filter_xml_array

subroutine read_xml_type_filter_xml( info, starttag, endtag, attribs, noattribs, data, nodata, &
              dvar, has_dvar )
   type(XML_PARSE)                                 :: info
   character(len=*), intent(in)                    :: starttag
   logical, intent(inout)                          :: endtag
   character(len=*), dimension(:,:), intent(inout) :: attribs
   integer, intent(inout)                          :: noattribs
   character(len=*), dimension(:), intent(inout)   :: data
   integer, intent(inout)                          :: nodata
   type(filter_xml), intent(inout)  :: dvar
   logical, intent(inout)                       :: has_dvar

   integer                                      :: att_
   integer                                      :: noatt_
   logical                                      :: error
   logical                                      :: endtag_org
   character(len=len(starttag))                 :: tag
   logical                                         :: has_type
   logical                                         :: has_bins
   logical                                         :: has_groups
   has_type                             = .false.
   has_bins                             = .false.
   has_groups                           = .false.
   call init_xml_type_filter_xml(dvar)
   has_dvar = .true.
   error  = .false.
   att_   = 0
   noatt_ = noattribs+1
   endtag_org = endtag
   do
      if ( nodata /= 0 ) then
         noattribs = 0
         tag = starttag
      elseif ( att_ < noatt_ .and. noatt_ > 1 ) then
         att_      = att_ + 1
         if ( att_ <= noatt_-1 ) then
            tag       = attribs(1,att_)
            data(1)   = attribs(2,att_)
            noattribs = 0
            nodata    = 1
            endtag    = .false.
         else
            tag       = starttag
            noattribs = 0
            nodata    = 0
            endtag    = .true.
            cycle
         endif
      else
         if ( endtag_org ) then
            return
         else
            call xml_get( info, tag, endtag, attribs, noattribs, data, nodata )
            if ( xml_error(info) ) then
               write(lurep_,*) 'Error reading input file!'
               error = .true.
               return
            endif
         endif
      endif
      if ( endtag .and. tag == starttag ) then
         exit
      endif
      if ( endtag .and. noattribs == 0 ) then
         if ( xml_ok(info) ) then
            cycle
         else
            exit
         endif
      endif
      select case( tag )
      case('type')
         call read_xml_word( &
            info, tag, endtag, attribs, noattribs, data, nodata, &
            dvar%type, has_type )
      case('bins')
         call read_xml_word_array( &
            info, tag, endtag, attribs, noattribs, data, nodata, &
            dvar%bins, has_bins )
      case('groups')
         call read_xml_word( &
            info, tag, endtag, attribs, noattribs, data, nodata, &
            dvar%groups, has_groups )
      case ('comment', '!--')
         ! Simply ignore
      case default
         if ( strict_ ) then
            error = .true.
            call xml_report_errors( info, &
               'Unknown or wrongly placed tag: ' // trim(tag))
         endif
      end select
      nodata = 0
      if ( .not. xml_ok(info) ) exit
   end do
   if ( .not. has_bins ) then
      has_dvar = .false.
      call xml_report_errors(info, 'Missing data on bins')
   endif
end subroutine read_xml_type_filter_xml
subroutine init_xml_type_filter_xml_array( dvar )
   type(filter_xml), dimension(:), pointer :: dvar
   if ( associated( dvar ) ) then
      deallocate( dvar )
   endif
   allocate( dvar(0) )
end subroutine init_xml_type_filter_xml_array
subroutine init_xml_type_filter_xml(dvar)
   type(filter_xml) :: dvar
   dvar%type = ''
   dvar%groups = ''
end subroutine init_xml_type_filter_xml
subroutine write_xml_type_filter_xml_array( &
      info, tag, indent, dvar )
   type(XML_PARSE)                                 :: info
   character(len=*), intent(in)                    :: tag
   integer                                         :: indent
   type(filter_xml), dimension(:)        :: dvar
   integer                                         :: i
   do i = 1,size(dvar)
       call write_xml_type_filter_xml( info, tag, indent, dvar(i) )
   enddo
end subroutine write_xml_type_filter_xml_array

subroutine write_xml_type_filter_xml( &
      info, tag, indent, dvar )
   type(XML_PARSE)                                 :: info
   character(len=*), intent(in)                    :: tag
   integer                                         :: indent
   type(filter_xml)                      :: dvar
   character(len=100)                              :: indentation
   indentation = ' '
   write(info%lun, '(4a)' ) indentation(1:min(indent,100)),&
       '<',trim(tag), '>'
   call write_to_xml_word( info, 'type', indent+3, dvar%type)
   call write_to_xml_word_array( info, 'bins', indent+3, dvar%bins)
   call write_to_xml_word( info, 'groups', indent+3, dvar%groups)
   write(info%lun,'(4a)') indentation(1:min(indent,100)), &
       '</' //trim(tag) // '>'
end subroutine write_xml_type_filter_xml

subroutine read_xml_type_tally_xml_array( &
      info, tag, endtag, attribs, noattribs, data, nodata, &
      dvar, has_dvar )
   type(XML_PARSE)                                 :: info
   character(len=*), intent(inout)                 :: tag
   logical, intent(inout)                          :: endtag
   character(len=*), dimension(:,:), intent(inout) :: attribs
   integer, intent(inout)                          :: noattribs
   character(len=*), dimension(:), intent(inout)   :: data
   integer, intent(inout)                          :: nodata
   type(tally_xml), dimension(:), pointer :: dvar
   logical, intent(inout)                       :: has_dvar

   integer                                      :: newsize
   type(tally_xml), dimension(:), pointer :: newvar

   newsize = size(dvar) + 1
   allocate( newvar(1:newsize) )
   newvar(1:newsize-1) = dvar
   deallocate( dvar )
   dvar => newvar

   call read_xml_type_tally_xml( info, tag, endtag, attribs, noattribs, data, nodata, &
              dvar(newsize), has_dvar )
end subroutine read_xml_type_tally_xml_array

subroutine read_xml_type_tally_xml( info, starttag, endtag, attribs, noattribs, data, nodata, &
              dvar, has_dvar )
   type(XML_PARSE)                                 :: info
   character(len=*), intent(in)                    :: starttag
   logical, intent(inout)                          :: endtag
   character(len=*), dimension(:,:), intent(inout) :: attribs
   integer, intent(inout)                          :: noattribs
   character(len=*), dimension(:), intent(inout)   :: data
   integer, intent(inout)                          :: nodata
   type(tally_xml), intent(inout)  :: dvar
   logical, intent(inout)                       :: has_dvar

   integer                                      :: att_
   integer                                      :: noatt_
   logical                                      :: error
   logical                                      :: endtag_org
   character(len=len(starttag))                 :: tag
   logical                                         :: has_id
   logical                                         :: has_label
   logical                                         :: has_estimator
   logical                                         :: has_filter
   logical                                         :: has_scores
   logical                                         :: has_nuclides
   logical                                         :: has_filters
   has_id                               = .false.
   has_label                            = .false.
   has_estimator                        = .false.
   has_filter                           = .false.
   allocate(dvar%filter(0))
   has_scores                           = .false.
   has_nuclides                         = .false.
   has_filters                          = .false.
   allocate(dvar%filters(0))
   call init_xml_type_tally_xml(dvar)
   has_dvar = .true.
   error  = .false.
   att_   = 0
   noatt_ = noattribs+1
   endtag_org = endtag
   do
      if ( nodata /= 0 ) then
         noattribs = 0
         tag = starttag
      elseif ( att_ < noatt_ .and. noatt_ > 1 ) then
         att_      = att_ + 1
         if ( att_ <= noatt_-1 ) then
            tag       = attribs(1,att_)
            data(1)   = attribs(2,att_)
            noattribs = 0
            nodata    = 1
            endtag    = .false.
         else
            tag       = starttag
            noattribs = 0
            nodata    = 0
            endtag    = .true.
            cycle
         endif
      else
         if ( endtag_org ) then
            return
         else
            call xml_get( info, tag, endtag, attribs, noattribs, data, nodata )
            if ( xml_error(info) ) then
               write(lurep_,*) 'Error reading input file!'
               error = .true.
               return
            endif
         endif
      endif
      if ( endtag .and. tag == starttag ) then
         exit
      endif
      if ( endtag .and. noattribs == 0 ) then
         if ( xml_ok(info) ) then
            cycle
         else
            exit
         endif
      endif
      select case( tag )
      case('id')
         call read_xml_integer( &
            info, tag, endtag, attribs, noattribs, data, nodata, &
            dvar%id, has_id )
      case('label')
         call read_xml_word( &
            info, tag, endtag, attribs, noattribs, data, nodata, &
            dvar%label, has_label )
      case('estimator')
         call read_xml_word( &
            info, tag, endtag, attribs, noattribs, data, nodata, &
            dvar%estimator, has_estimator )
      case('filter')
         call read_xml_type_filter_xml_array( &
            info, tag, endtag, attribs, noattribs, data, nodata, &
            dvar%filter, has_filter )
      case('scores')
         call read_xml_word_array( &
            info, tag, endtag, attribs, noattribs, data, nodata, &
            dvar%scores, has_scores )
      case('nuclides')
         call read_xml_word_array( &
            info, tag, endtag, attribs, noattribs, data, nodata, &
            dvar%nuclides, has_nuclides )
      case('filters')
         call read_xml_type_filter_xml_array( &
            info, tag, endtag, attribs, noattribs, data, nodata, &
            dvar%filters, has_filters )
      case ('comment', '!--')
         ! Simply ignore
      case default
         if ( strict_ ) then
            error = .true.
            call xml_report_errors( info, &
               'Unknown or wrongly placed tag: ' // trim(tag))
         endif
      end select
      nodata = 0
      if ( .not. xml_ok(info) ) exit
   end do
   if ( .not. has_id ) then
      has_dvar = .false.
      call xml_report_errors(info, 'Missing data on id')
   endif
   if ( .not. has_filter ) then
      has_dvar = .false.
      call xml_report_errors(info, 'Missing data on filter')
   endif
   if ( .not. has_scores ) then
      has_dvar = .false.
      call xml_report_errors(info, 'Missing data on scores')
   endif
   if ( .not. has_nuclides ) then
      has_dvar = .false.
      call xml_report_errors(info, 'Missing data on nuclides')
   endif
   if ( .not. has_filters ) then
      has_dvar = .false.
      call xml_report_errors(info, 'Missing data on filters')
   endif
end subroutine read_xml_type_tally_xml
subroutine init_xml_type_tally_xml_array( dvar )
   type(tally_xml), dimension(:), pointer :: dvar
   if ( associated( dvar ) ) then
      deallocate( dvar )
   endif
   allocate( dvar(0) )
end subroutine init_xml_type_tally_xml_array
subroutine init_xml_type_tally_xml(dvar)
   type(tally_xml) :: dvar
   dvar%label = ''
   dvar%estimator = ''
end subroutine init_xml_type_tally_xml
subroutine write_xml_type_tally_xml_array( &
      info, tag, indent, dvar )
   type(XML_PARSE)                                 :: info
   character(len=*), intent(in)                    :: tag
   integer                                         :: indent
   type(tally_xml), dimension(:)        :: dvar
   integer                                         :: i
   do i = 1,size(dvar)
       call write_xml_type_tally_xml( info, tag, indent, dvar(i) )
   enddo
end subroutine write_xml_type_tally_xml_array

subroutine write_xml_type_tally_xml( &
      info, tag, indent, dvar )
   type(XML_PARSE)                                 :: info
   character(len=*), intent(in)                    :: tag
   integer                                         :: indent
   type(tally_xml)                      :: dvar
   character(len=100)                              :: indentation
   indentation = ' '
   write(info%lun, '(4a)' ) indentation(1:min(indent,100)),&
       '<',trim(tag), '>'
   call write_to_xml_integer( info, 'id', indent+3, dvar%id)
   call write_to_xml_word( info, 'label', indent+3, dvar%label)
   call write_to_xml_word( info, 'estimator', indent+3, dvar%estimator)
   call write_xml_type_filter_xml_array( info, 'filter', indent+3, dvar%filter)
   call write_to_xml_word_array( info, 'scores', indent+3, dvar%scores)
   call write_to_xml_word_array( info, 'nuclides', indent+3, dvar%nuclides)
   call write_xml_type_filter_xml_array( info, 'filters', indent+3, dvar%filters)
   write(info%lun,'(4a)') indentation(1:min(indent,100)), &
       '</' //trim(tag) // '>'
end subroutine write_xml_type_tally_xml

subroutine read_xml_file_tallies_t(fname, lurep, errout)
   character(len=*), intent(in)           :: fname
   integer, intent(in), optional          :: lurep
   logical, intent(out), optional         :: errout

   type(XML_PARSE)                        :: info
   logical                                :: error
   character(len=80)                      :: tag
   character(len=80)                      :: starttag
   logical                                :: endtag
   character(len=250), dimension(1:2,1:20) :: attribs
   integer                                :: noattribs
   character(len=1000), dimension(1:1000)  :: data
   integer                                :: nodata
   logical                                         :: has_mesh_
   logical                                         :: has_tally_
   logical                                         :: has_separate_
   has_mesh_                            = .false.
   allocate(mesh_(0))
   has_tally_                           = .false.
   allocate(tally_(0))
   has_separate_                        = .false.

   call init_xml_file_tallies_t
   call xml_open( info, fname, .true. )
   call xml_options( info, report_errors=.false., ignore_whitespace=.true.)
   lurep_ = 0
   if ( present(lurep) ) then
      lurep_ = lurep
      call xml_options( info, report_lun=lurep )
   endif
   do
      call xml_get( info, starttag, endtag, attribs, noattribs, &
         data, nodata)
      if ( starttag /= '!--' ) exit
   enddo
   if ( starttag /= "tallies" ) then
      call xml_report_errors( info, &
         'XML-file should have root element "tallies"')
      error = .true.
      call xml_close(info)
      return
   endif
   strict_ = .false.
   error = .false.
   do
      call xml_get( info, tag, endtag, attribs, noattribs, data, nodata )
      if ( xml_error(info) ) then
         write(lurep_,*) 'Error reading input file!'
         error = .true.
         return
      endif
      if ( endtag .and. tag == starttag ) then
         exit
      endif
      if ( endtag .and. noattribs == 0 ) then
         if ( xml_ok(info) ) then
            cycle
         else
            exit
         endif
      endif
      select case( tag )
      case('mesh')
         call read_xml_type_mesh_xml_array( &
            info, tag, endtag, attribs, noattribs, data, nodata, &
            mesh_, has_mesh_ )
      case('tally')
         call read_xml_type_tally_xml_array( &
            info, tag, endtag, attribs, noattribs, data, nodata, &
            tally_, has_tally_ )
      case('assume_separate')
         call read_xml_word( &
            info, tag, endtag, attribs, noattribs, data, nodata, &
            separate_, has_separate_ )
      case ('comment', '!--')
         ! Simply ignore
      case default
         if ( strict_ ) then
            error = .true.
            call xml_report_errors( info, &
               'Unknown or wrongly placed tag: ' // trim(tag))
         endif
      end select
      nodata = 0
      if ( .not. xml_ok(info) ) exit
   end do
   if ( .not. has_mesh_ ) then
      error = .true.
      call xml_report_errors(info, 'Missing data on mesh_')
   endif
   if ( .not. has_tally_ ) then
      error = .true.
      call xml_report_errors(info, 'Missing data on tally_')
   endif
   if ( .not. has_separate_ ) then
      error = .true.
      call xml_report_errors(info, 'Missing data on separate_')
   endif
   if ( present(errout) ) errout = error
   call xml_close(info)
end subroutine

subroutine write_xml_file_tallies_t(fname, lurep)
   character(len=*), intent(in)           :: fname
   integer, intent(in), optional          :: lurep

   type(XML_PARSE)                        :: info
   integer                                :: indent = 0

   call xml_open( info, fname, .false. )
   call xml_options( info, report_errors=.true.)
   if ( present(lurep) ) then
       call xml_options( info, report_errors=.true.)
   endif
   write(info%lun,'(a)') &
      '<tallies>'
   call write_xml_type_mesh_xml_array( info, 'mesh', indent+3, mesh_)
   call write_xml_type_tally_xml_array( info, 'tally', indent+3, tally_)
   call write_to_xml_word( info, 'assume_separate', indent+3, separate_)
   write(info%lun,'(a)') '</tallies>'
   call xml_close(info)
end subroutine

subroutine init_xml_file_tallies_t

end subroutine

end module
